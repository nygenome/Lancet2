#ifndef SRC_LANCET_CBDG_GRAPH_H_
#define SRC_LANCET_CBDG_GRAPH_H_

#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/graph_complexity.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/read.h"
#include "lancet/cbdg/traversal_index.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lancet::cbdg {

class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeTable = absl::flat_hash_map<NodeID, NodePtr>;
  using RegionPtr = std::shared_ptr<hts::Reference::Region const>;
  using ReadList = absl::Span<Read const>;

  static constexpr usize DEFAULT_MIN_KMER_LEN = 13;
  static constexpr usize DEFAULT_MAX_KMER_LEN = 127;
  static constexpr usize MAX_ALLOWED_KMER_LEN = 255;

  static constexpr u32 DEFAULT_MIN_NODE_COV = 2;
  static constexpr u32 DEFAULT_MIN_ANCHOR_COV = 5;
  static constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1'048'576;

  static constexpr u16 DEFAULT_KMER_STEP_LEN = 6;

  struct Params {
    std::filesystem::path mOutGraphsDir;

    usize mMinKmerLen = DEFAULT_MIN_KMER_LEN;
    usize mMaxKmerLen = DEFAULT_MAX_KMER_LEN;

    u32 mMinNodeCov = DEFAULT_MIN_NODE_COV;
    u32 mMinAnchorCov = DEFAULT_MIN_ANCHOR_COV;

    u16 mKmerStepLen = DEFAULT_KMER_STEP_LEN;
  };

  explicit Graph(Params params) : mParams(std::move(params)) {}

  [[nodiscard]] auto CurrentK() const noexcept -> usize { return mCurrK; }

  /// Const access to the node table — used by graph_complexity.cpp free functions.
  [[nodiscard]] auto Nodes() const noexcept -> NodeTable const& { return mNodes; }

  class Path {
   public:
    Path() = default;

    void AppendSequence(std::string_view seq);
    void ReserveSequence(usize len) { mSequence.reserve(len); }
    void AddNodeCoverage(u32 cov);
    void Finalize();

    [[nodiscard]] auto IsEmpty() const -> bool { return mSequence.empty(); }
    [[nodiscard]] auto Sequence() const -> std::string_view { return mSequence; }
    /// Average read coverage across all nodes constituting this path
    [[nodiscard]] auto MeanCoverage() const -> f64 { return mMeanCov; }
    /// Median read coverage across all nodes constituting this path
    [[nodiscard]] auto MedianCoverage() const -> f64 { return mMedianCov; }
    /// Standard deviation of read coverage across nodes
    [[nodiscard]] auto StdDevCoverage() const -> f64 { return mStdDevCov; }
    /// Coefficient of Variation (StdDev / Mean) -> higher means more coverage fluctuation
    [[nodiscard]] auto CoefficientOfVariationCoverage() const -> f64 { return mCvCov; }
    /// Quartile Coefficient of Variation ((Q3 - Q1) / (Q3 + Q1)) -> robust against outliers
    [[nodiscard]] auto QuartileCoefficientOfVariation() const -> f64 { return mQCvCov; }
    /// Aggregate total of all node coverages on this path
    [[nodiscard]] auto TotalCoverage() const -> f64 { return mTotalCov; }

   private:
    std::string mSequence;
    absl::InlinedVector<u32, 256> mNodeCoverages;
    f64 mMeanCov = 0.0;
    f64 mMedianCov = 0.0;
    f64 mStdDevCov = 0.0;
    f64 mCvCov = 0.0;
    f64 mQCvCov = 0.0;
    f64 mTotalCov = 0.0;
  };

  // First is always ref hap. Subsequent ALT haplotypes are sorted by descending
  // MeanCoverage, establishing Greedy Insertion Bias in downstream SPOA MSA.
  using CompHaps = std::vector<Path>;
  using GraphHaps = std::vector<CompHaps>;

  struct Result {
    GraphHaps mGraphHaplotypes;
    std::vector<usize> mAnchorStartIdxs;
    std::vector<GraphComplexity> mComponentMetrics;
  };

  [[nodiscard]] auto BuildComponentHaplotypes(RegionPtr region, ReadList reads) -> Result;

 private:
  usize mCurrK = 0;
  RegionPtr mRegion;
  ReadList mReads;
  NodeTable mNodes;
  Params mParams;

  std::vector<NodeID> mRefNodeIds;
  NodeIDPair mSourceAndSinkIds = {0, 0};

  using EdgeSet = absl::flat_hash_set<Edge>;
  using NodeIdSet = absl::flat_hash_set<NodeID>;

  void CompressGraph(usize component_id);
  void CompressNode(NodeID nid, Kmer::Ordering ord, NodeIdSet& compressed_ids) const;

  [[nodiscard]] auto FindCompressibleEdge(Node const& src, Kmer::Ordering ord) const
      -> std::optional<Edge>;
  [[nodiscard]] auto IsPotentialBuddyEdge(Node const& src, Edge const& conn) const -> bool;

  void RemoveTips(usize component_id);

  struct RefAnchor {
    // NOLINTNEXTLINE(cppcoreguidelines-use-enum-class)
    enum Kind : bool { SOURCE = true, SINK = false };

    NodeID mAnchorId = 0;
    usize mRefOffset = 0;
    bool mFoundAnchor = false;
  };

  [[nodiscard]] auto FindSource(usize component_id) const -> RefAnchor;
  [[nodiscard]] auto FindSink(usize component_id) const -> RefAnchor;

  /// O(V+E) cycle detection using three-color DFS on the flat adjacency list.
  /// Tracks (NodeIdx, Sign) state pairs — a node reached via '+' and via '-' are
  /// different states, correctly handling bidirected sign continuity (BCALM2).
  /// Returns true if any cycle is reachable from the source node.
  [[nodiscard]] static auto HasCycle(TraversalIndex const& idx) -> bool;

  /// Build a flat adjacency list from the frozen (fully-pruned) graph for a single
  /// component. The index maps NodeID -> contiguous u32, enabling O(1) array-based
  /// traversal state tracking. Built once, consumed by HasCycle and MaxFlow.
  [[nodiscard]] auto BuildTraversalIndex(usize component_id) const -> TraversalIndex;

  // ============================================================================
  // Inner Loop Modularity Methods
  // ============================================================================
  void PruneComponent(usize component_id);
  [[nodiscard]] auto EnumerateAndSortHaplotypes(usize comp_id, TraversalIndex const& trav_idx,
                                                std::string_view ref_anchor_seq) const
      -> std::vector<Path>;

  [[nodiscard]] auto ComputeComponentComplexity(usize component_id) const -> GraphComplexity;

  struct ComponentInfo {
    f64 mPctNodes = 0.0;
    usize mCompId = 0;
    usize mNumNodes = 0;
  };
  [[nodiscard]] auto MarkConnectedComponents() -> std::vector<ComponentInfo>;

  void RemoveLowCovNodes(usize component_id);
  void RemoveNode(NodeTable::iterator itr);
  void RemoveNodes(absl::Span<NodeID const> node_ids);

  // mateMer -> readName + sampleLabel, kmerHash
  struct MateMer {
    std::string_view mQname;
    Label::Tag mTagKind;
    u64 mKmerHash;

    auto operator==(MateMer const& rhs) const -> bool {
      return mQname == rhs.mQname && mTagKind == rhs.mTagKind && mKmerHash == rhs.mKmerHash;
    }

    template <typename H>
    friend auto AbslHashValue(H state, MateMer const& mmer) -> H {
      return H::combine(std::move(state), mmer.mQname, mmer.mTagKind, mmer.mKmerHash);
    }
  };
  void BuildGraph(absl::flat_hash_set<MateMer>& mate_mers);
  auto AddNodes(std::string_view sequence, Label label) -> std::vector<Node*>;

  [[nodiscard]] static auto HasExactOrApproxRepeat(std::string_view seq, usize window) -> bool;
  [[nodiscard]] static auto RefAnchorLength(RefAnchor const& source, RefAnchor const& sink,
                                            usize currk) -> usize;

  // NOLINTNEXTLINE(cppcoreguidelines-use-enum-class)
  enum State : u8 {
    FIRST_LOW_COV_REMOVAL = 0,
    FOUND_REF_ANCHORS = 1,
    FIRST_COMPRESSION = 2,
    SECOND_LOW_COV_REMOVAL = 3,
    SECOND_COMPRESSION = 4,
    SHORT_TIP_REMOVAL = 5,
    SHORT_LINK_REMOVAL = 6,
    FULLY_PRUNED_GRAPH = 7
  };

  void WriteDot([[maybe_unused]] State state, usize comp_id);

#ifdef LANCET_DEVELOP_MODE
  [[nodiscard]] static auto ToString(State state) -> std::string;

  template <class... Args>
  inline constexpr void WriteDotDevelop(Args&&... args) {
    WriteDot(std::forward<Args>(args)...);
  }
#else
  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  constexpr void WriteDotDevelop([[maybe_unused]] Args&&... /*unused*/) {}
#endif

  static void SerializeToDot(NodeTable const& graph, std::filesystem::path const& out_path,
                             usize comp_id = 0, NodeIdSet const& nodes_highlight = {},
                             EdgeSet const& edges_highlight = {},
                             NodeIdSet const& nodes_background = {},
                             EdgeSet const& edges_background = {});
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_H_

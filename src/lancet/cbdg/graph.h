#ifndef SRC_LANCET_CBDG_GRAPH_H_
#define SRC_LANCET_CBDG_GRAPH_H_

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/cbdg/edge.h"
#include "lancet/cbdg/kmer.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/node.h"
#include "lancet/cbdg/read.h"
#include "lancet/cbdg/traversal_index.h"
#include "lancet/hts/reference.h"

namespace lancet::cbdg {

class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeTable = absl::flat_hash_map<NodeID, NodePtr>;
  using RegionPtr = std::shared_ptr<const hts::Reference::Region>;
  using ReadList = absl::Span<const Read>;

  static constexpr usize DEFAULT_MIN_KMER_LEN = 13;
  static constexpr usize DEFAULT_MAX_KMER_LEN = 127;
  static constexpr usize MAX_ALLOWED_KMER_LEN = 255;

  static constexpr u32 DEFAULT_MIN_NODE_COV = 2;
  static constexpr u32 DEFAULT_MIN_ANCHOR_COV = 5;
  static constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1048576;

  static constexpr u16 DEFAULT_KMER_STEP_LEN = 6;

  struct Params {
    std::filesystem::path mOutGraphsDir;

    usize mMinKmerLen = DEFAULT_MIN_KMER_LEN;
    usize mMaxKmerLen = DEFAULT_MAX_KMER_LEN;

    u32 mMinNodeCov = DEFAULT_MIN_NODE_COV;
    u32 mMinAnchorCov = DEFAULT_MIN_ANCHOR_COV;

    u16 mKmerStepLen = DEFAULT_KMER_STEP_LEN;
  };

  Graph(Params params) : mParams(std::move(params)) {}

  [[nodiscard]] auto CurrentK() const noexcept -> usize { return mCurrK; }

  // First is always ref hap, rest are alts
  using CompHaps = std::vector<std::string>;
  using GraphHaps = std::vector<CompHaps>;

  struct Result {
    GraphHaps mGraphHaplotypes;
    std::vector<usize> mAnchorStartIdxs;
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

  [[nodiscard]] auto FindCompressibleEdge(const Node& src, Kmer::Ordering ord) const -> std::optional<Edge>;
  [[nodiscard]] auto IsPotentialBuddyEdge(const Node& src, const Edge& conn) const -> bool;

  void RemoveTips(usize component_id);

  struct RefAnchor {
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
  [[nodiscard]] auto HasCycle(const TraversalIndex& idx) const -> bool;

  /// Build a flat adjacency list from the frozen (fully-pruned) graph for a single
  /// component. The index maps NodeID -> contiguous u32, enabling O(1) array-based
  /// traversal state tracking. Built once, consumed by HasCycle and MaxFlow.
  [[nodiscard]] auto BuildTraversalIndex(usize component_id) const -> TraversalIndex;

  /// Lightweight O(V+E) graph complexity metrics for debug logging and
  /// pathological graph detection. Computed right before MaxFlow.
  struct GraphComplexity {
    usize mNumNodes = 0;                ///< V (nodes in this component)
    usize mNumEdges = 0;                ///< E (unduplicated edge count)
    usize mCyclomaticComplexity = 0;    ///< M = E - V + 1 (single component)
    f64 mEdgeToNodeDensity = 0.0;       ///< E / V
    usize mMaxSingleDirDegree = 0;      ///< max outgoing edges in any single sign direction
    usize mNumBranchPoints = 0;         ///< nodes with >= 2 edges in some direction
    /// Thresholds for classifying a graph as too complex for efficient walk
    /// enumeration at the current k-mer size. When IsComplex() returns true,
    /// the caller retries with a larger k to collapse branches.
    ///
    /// Biological interpretation (in a ~1000bp window with k=25, ~975 nodes):
    ///   CC = 50 ≈ 50 independent variant signals in 1000bp (~1 per 20bp),
    ///   which is ~50× the normal human heterozygosity rate (~1 SNP/1000bp).
    ///   Br = 50 ≈ 5% of all positions forking into multiple continuations.
    ///   Together, CC>=50 AND Br>=50 indicates the graph is not driven by
    ///   clean biological variation but by tandem repeats (STRs/VNTRs),
    ///   segmental duplications, structural variant breakpoints, or mapping
    ///   artifacts. A larger k collapses short repeats and merges paralogous
    ///   paths, simplifying the graph toward its true biological structure.
    ///
    /// Derived from whole-chromosome profiling (chr4, 233,479 windows, 235,533
    /// component entries).
    ///
    ///   Threshold analysis (AND — both conditions must hold):
    ///     CC>=50  AND Br>=50:   146 wins, avg  5.8s | 233,333 normal, avg 414ms (14×)
    ///     CC>=50  AND Br>=100:   76 wins, avg 10.4s | 233,403 normal, avg 414ms (25×)
    ///     CC>=100 AND Br>=100:   31 wins, avg 15.2s | 233,448 normal, avg 415ms (37×)
    ///     CC>=100 AND Br>=200:   18 wins, avg 25.0s | 233,461 normal, avg 415ms (60×)
    ///
    ///   Using AND (vs OR) avoids skipping windows where only one metric is
    ///   elevated — e.g. a high-CC but low-branch graph may still enumerate
    ///   quickly. Both metrics must agree that the graph is pathological.
    ///
    ///   Distribution stats (235,533 components):
    ///     CC:  median=3, p95=9, p99=16, max=487
    ///     Br:  median=6, p95=32, p99=32, max=1127
    ///
    static constexpr usize MAX_CYCLOMATIC_COMPLEXITY = 50;
    static constexpr usize MAX_BRANCH_POINTS = 50;

    /// Returns true if the graph is too complex for efficient walk enumeration.
    /// When true, the caller should retry with a larger k-mer size to collapse
    /// branches and reduce cyclomatic complexity.
    [[nodiscard]] constexpr auto IsComplex() const -> bool {
      return mCyclomaticComplexity >= MAX_CYCLOMATIC_COMPLEXITY && mNumBranchPoints >= MAX_BRANCH_POINTS;
    }
  };
  [[nodiscard]] auto ComputeGraphComplexity(usize component_id) const -> GraphComplexity;

  struct ComponentInfo {
    f64 mPctNodes = 0.0;
    usize mCompId = 0;
    usize mNumNodes = 0;
  };
  [[nodiscard]] auto MarkConnectedComponents() -> std::vector<ComponentInfo>;

  void RemoveLowCovNodes(usize component_id);
  void RemoveNode(NodeTable::iterator itr);
  void RemoveNodes(absl::Span<const NodeID> node_ids);

  // mateMer -> readName + sampleLabel, kmerHash
  struct MateMer {
    std::string_view mQname;
    Label::Tag mTagKind;
    u64 mKmerHash;

    auto operator==(const MateMer& rhs) const -> bool {
      return mQname == rhs.mQname && mTagKind == rhs.mTagKind && mKmerHash == rhs.mKmerHash;
    }

    template <typename H>
    friend auto AbslHashValue(H state, const MateMer& mm) -> H {
      return H::combine(std::move(state), mm.mQname, mm.mTagKind, mm.mKmerHash);
    }
  };
  void BuildGraph(absl::flat_hash_set<MateMer>& mate_mers);
  auto AddNodes(std::string_view sequence, Label label) -> std::vector<Node*>;

  [[nodiscard]] static auto HasExactOrApproxRepeat(std::string_view seq, usize window) -> bool;
  [[nodiscard]] static auto RefAnchorLength(const RefAnchor& source, const RefAnchor& sink, usize currk) -> usize;

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
  constexpr inline void WriteDotDevelop(Args&&... args) {
    WriteDot(std::forward<Args>(args)...);
  }
#else
  template <class... Args>
  // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
  constexpr void WriteDotDevelop([[maybe_unused]] Args&&... /*unused*/) {}
#endif

  static void SerializeToDot(const NodeTable& graph, const std::filesystem::path& out_path, usize comp_id = 0,
                             const NodeIdSet& nodes_highlight = {}, const EdgeSet& edges_highlight = {},
                             const NodeIdSet& nodes_background = {}, const EdgeSet& edges_background = {});
};

}  // namespace lancet::cbdg

#endif  // SRC_LANCET_CBDG_GRAPH_H_

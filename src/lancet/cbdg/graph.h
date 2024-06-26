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
#include "lancet/hts/reference.h"

namespace lancet::cbdg {

class Graph {
 public:
  using NodePtr = std::unique_ptr<Node>;
  using NodeTable = absl::flat_hash_map<NodeID, NodePtr>;
  using RegionPtr = std::shared_ptr<const hts::Reference::Region>;
  using ReadList = absl::Span<const Read>;

  static constexpr usize DEFAULT_MIN_KMER_LEN = 31;
  static constexpr usize DEFAULT_MAX_KMER_LEN = 133;
  static constexpr usize MAX_ALLOWED_KMER_LEN = 255;

  static constexpr u32 DEFAULT_MIN_NODE_COV = 2;
  static constexpr u32 DEFAULT_MIN_ANCHOR_COV = 5;
  static constexpr u32 DEFAULT_GRAPH_TRAVERSAL_LIMIT = 1e6;

  static constexpr u16 DEFAULT_KMER_STEP_LEN = 4;

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

  [[nodiscard]] auto HasCycle() const -> bool;
  void HasCycle(const Node& node, NodeIdSet& traversed, bool& found_cycle, usize& recursion_depth) const;

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
  using MateMer = std::pair<std::string, u64>;
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

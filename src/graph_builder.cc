#include "lancet2/graph_builder.h"

#include <utility>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#include "lancet2/assert_macro.h"
#include "lancet2/canonical_kmers.h"
#include "lancet2/core_enums.h"
#include "lancet2/kmer.h"
#include "lancet2/log_macros.h"
#include "lancet2/utils.h"
#include "spdlog/spdlog.h"

#ifndef NDEBUG
#include "lancet2/timer.h"
#endif

namespace lancet2 {
GraphBuilder::GraphBuilder(std::shared_ptr<const RefWindow> w, absl::Span<const ReadInfo> reads, double avg_cov,
                           std::shared_ptr<const CliParams> p)
    : avgCov(avg_cov), window(std::move(w)), params(std::move(p)), sampleReads(reads) {}

auto GraphBuilder::BuildGraph(usize min_k, usize max_k) -> std::unique_ptr<Graph> {
#ifndef NDEBUG
  Timer timer;
#endif

  const auto windowId = window->ToRegionString();
  LOG_DEBUG("Starting to build graph for {} using minK={}", windowId, min_k);

  for (currentK = min_k; currentK <= max_k; currentK += 2) {
    if (utils::HasRepeatKmer(window->GetSeqView(), currentK)) continue;
    if (utils::HasAlmostRepeatKmer(window->GetSeqView(), currentK, params->maxRptMismatch)) continue;

    nodesMap.clear();
    BuildSampleNodes();
    BuildRefNodes();
    BuildNode(MOCK_SOURCE_ID);
    BuildNode(MOCK_SINK_ID);

    if (params->kmerRecoveryOn) RecoverKmers();
    break;
  }

#ifndef NDEBUG
  LOG_DEBUG("Built graph for {} with K={} | Runtime={}", windowId, currentK, timer.HumanRuntime());
#endif
  return std::make_unique<Graph>(window, std::move(nodesMap), avgCov, currentK, params);
}

void GraphBuilder::BuildSampleNodes() {
  // mateMer -> readName + sampleLabel, kmerHash
  using MateMer = std::pair<std::string, u64>;
  absl::flat_hash_set<MateMer> seenMateMers;
  LANCET_ASSERT(!sampleReads.empty());  // NOLINT

  for (const auto& rd : sampleReads) {
    const auto result = BuildNodes(rd.sequence);
    const auto qualMers = KMovingSubstrs(rd.quality, currentK);
    LANCET_ASSERT(result.nodeIDs.size() == qualMers.size());  // NOLINT

    for (usize idx = 0; idx < result.nodeIDs.size() - 1; idx++) {
      const auto firstId = result.nodeIDs[idx];
      const auto secondId = result.nodeIDs[idx + 1];

      auto itr1 = nodesMap.find(firstId);
      auto itr2 = nodesMap.find(secondId);
      LANCET_ASSERT(itr1 != nodesMap.end());  // NOLINT
      LANCET_ASSERT(itr2 != nodesMap.end());  // NOLINT

      const auto firstEk = MakeEdgeKind(itr1->second->GetOrientation(), itr2->second->GetOrientation());
      itr1->second->EmplaceEdge(secondId, firstEk);
      itr2->second->EmplaceEdge(firstId, ReverseEdgeKind(firstEk));

      const auto nodeLabel = rd.label == SampleLabel::TUMOR ? KmerLabel::TUMOR : KmerLabel::NORMAL;
      itr1->second->UpdateLabel(nodeLabel);
      itr2->second->UpdateLabel(nodeLabel);

      // pick node to update. make sure count info is not updated twice for same node
      // since we are moving with consecutive kmers, we know unique nodes by iteration index
      auto currItr = idx == 0 ? itr1 : itr2;
      const auto& currQual = idx == 0 ? qualMers[idx] : qualMers[idx + 1];
      const auto isCurrFwd = currItr->second->GetOrientation() == Strand::FWD;

      currItr->second->UpdateQual(isCurrFwd ? currQual : std::string(currQual.crbegin(), currQual.crend()));
      if (params->tenxMode) currItr->second->UpdateHPInfo(rd, params->minBaseQual);

      // Add label to key, so that keys are unique for tumor and normal samples
      const auto mmId = std::make_pair(rd.readName + ToString(rd.label), currItr->first);
      const auto isUniqMateMer = seenMateMers.find(mmId) == seenMateMers.end();
      if (isUniqMateMer) {
        currItr->second->UpdateCovInfo(rd, params->minBaseQual, params->tenxMode);
        seenMateMers.insert(mmId);
      }
    }
  }

  const auto windowId = window->ToRegionString();
  LOG_DEBUG("Combined sample coverage for {} is ~{:.2f}x", windowId, avgCov);
  LOG_DEBUG("Built {} sample nodes in graph for {} using K={}", nodesMap.size(), windowId, currentK);
}

void GraphBuilder::BuildRefNodes() {
  const auto refDataBuilt = refNmlData.size() == params->windowLength && refTmrData.size() == params->windowLength;
  const auto refMerHashes = CanonicalKmerHashes(window->GetSeqView(), currentK);

  if (!refDataBuilt) BuildRefData(absl::MakeConstSpan(refMerHashes));
  usize numMarked = 0;

  for (usize idx = 0; idx < refMerHashes.size() - 1; idx++) {
    auto itr1 = nodesMap.find(refMerHashes[idx]);
    const auto foundNode1 = itr1 != nodesMap.end();

    auto itr2 = nodesMap.find(refMerHashes[idx + 1]);
    const auto foundNode2 = itr2 != nodesMap.end();

    if (foundNode1) {
      numMarked++;
      itr1->second->UpdateLabel(KmerLabel::REFERENCE);
    }

    if (foundNode2) {
      numMarked++;
      itr2->second->UpdateLabel(KmerLabel::REFERENCE);
    }

    if (foundNode1 && foundNode2) {
      const auto firstEk = MakeEdgeKind(itr1->second->GetOrientation(), itr2->second->GetOrientation());
      itr1->second->EmplaceEdge(itr2->first, firstEk);
      itr2->second->EmplaceEdge(itr1->first, ReverseEdgeKind(firstEk));
    }
  }

  LOG_DEBUG("Marked {} existing nodes as reference in graph for {} with K={}", numMarked, window->ToRegionString(),
            currentK);
}

auto GraphBuilder::BuildNodes(absl::string_view seq) -> GraphBuilder::BuildNodesResult {
  const auto mers = CanonicalKmers(seq, currentK);
  BuildNodesResult result{std::vector<NodeIdentifier>(), 0, mers.size()};
  result.nodeIDs.reserve(mers.size());

  for (auto&& mer : mers) {
    const auto merId = mer.GetHash();
    result.nodeIDs.emplace_back(merId);

    if (nodesMap.find(merId) == nodesMap.end()) {
      nodesMap.try_emplace(merId, std::make_unique<Node>(mer));
      result.numNodesBuilt++;
    }
  }

  return result;
}

auto GraphBuilder::BuildNode(NodeIdentifier node_id) -> GraphBuilder::BuildNodeResult {
  if (nodesMap.find(node_id) == nodesMap.end()) {
    nodesMap.try_emplace(node_id, std::make_unique<Node>(node_id));
    return BuildNodeResult{node_id, true};
  }

  return BuildNodeResult{node_id, false};
}

void GraphBuilder::RecoverKmers() {
  constexpr u16 minReadSupport = 2;
  usize numRecovered = 0;

  for (Graph::NodeContainer::const_reference p : nodesMap) {
    // recover only tumor singletons
    if (p.second->SampleCount(SampleLabel::TUMOR) != 1) continue;

    // identify positions with low quality bases
    const auto lowQualBits = p.second->LowQualPositions(params->minBaseQual);
    for (usize basePos = 0; basePos < lowQualBits.size(); ++basePos) {
      if (!lowQualBits[basePos]) continue;

      // get 6 alternative kmers within 1 edit distance to original kmer
      // TODO(omicsnut): pick alt kmer with highest/lowest support and only increment that?
      for (const auto& altKmer : MutateSeq(p.second->GetSeqView(), basePos)) {
        // check if alternative kmer exists in the graph
        auto altKmerItr = nodesMap.find(Kmer(altKmer).GetHash());
        if (altKmerItr == nodesMap.end()) continue;
        // check if node with alternative kmer has atleast one HQ tumor base at `basePos`
        // and has atleast `minReadSupport` tumor reads supporting it
        auto& altNode = altKmerItr->second;
        const auto hqCovAtPos = altNode->CovAt(SampleLabel::TUMOR, basePos).BQPassTotalCov();
        const auto readSupport = altNode->SampleCount(SampleLabel::TUMOR);
        if (hqCovAtPos <= 0 || readSupport < minReadSupport) continue;

        // TODO(omicsnut): why do we pick FWD when non 0, why not lower/higher?
        const auto tumorFwdCount = altNode->SampleCount(SampleLabel::TUMOR, Strand::FWD);
        const auto strandToIncrement = tumorFwdCount > 0 ? Strand::FWD : Strand::REV;
        altNode->IncrementCov(SampleLabel::TUMOR, strandToIncrement, basePos);
        numRecovered++;
      }
    }
  }

  if (numRecovered > 0) {
    LOG_DEBUG("Recovered %d singleton tumor %d-mers from graph for %s", numRecovered, currentK,
              window->ToRegionString());
  }
}

void GraphBuilder::BuildRefData(absl::Span<const usize> ref_mer_hashes) {
  NodeCov refCovs;
  NodeHP refHPs;

  refCovs.Reserve(params->windowLength);
  if (params->tenxMode) refHPs.Reserve(params->windowLength);

  for (const auto& merHash : ref_mer_hashes) {
    const auto itr = nodesMap.find(merHash);
    const auto foundNode = itr != nodesMap.end();
    const auto shouldReverse = foundNode && itr->second->GetOrientation() == Strand::REV;

    if (refCovs.IsEmpty()) {
      refCovs.Reserve(params->windowLength);
      if (params->tenxMode) refHPs.Reserve(params->windowLength);

      refCovs = foundNode ? itr->second->CovData() : NodeCov(currentK);
      if (params->tenxMode) refHPs = foundNode ? itr->second->HPData() : NodeHP(refCovs);
      continue;
    }

    const auto currCov = foundNode ? itr->second->CovData() : NodeCov(currentK);
    refCovs.MergeBuddy(currCov, BuddyPosition::FRONT, shouldReverse, currentK);

    if (params->tenxMode) {
      const auto currHp = foundNode ? itr->second->HPData() : NodeHP(currCov);
      refHPs.MergeBuddy(currHp, BuddyPosition::FRONT, shouldReverse, currentK);
    }
  }

  refNmlData = params->tenxMode
                   ? BuildHPCovs(refCovs.BaseCovs(SampleLabel::NORMAL), refHPs.BaseHPs(SampleLabel::NORMAL))
                   : BuildHPCovs(refCovs.BaseCovs(SampleLabel::NORMAL));

  refTmrData = params->tenxMode ? BuildHPCovs(refCovs.BaseCovs(SampleLabel::TUMOR), refHPs.BaseHPs(SampleLabel::TUMOR))
                                : BuildHPCovs(refCovs.BaseCovs(SampleLabel::TUMOR));

  LANCET_ASSERT(refNmlData.size() == params->windowLength + 1);  // NOLINT
  LANCET_ASSERT(refTmrData.size() == params->windowLength + 1);  // NOLINT
}

auto GraphBuilder::MutateSeq(absl::string_view seq, usize base_pos) -> std::vector<std::string> {
  std::vector<std::string> result;
  constexpr usize numMutatedSeqs = 6;
  result.reserve(numMutatedSeqs);  // 3 alt-kmers + 3 reverse complements
  static constexpr std::array<char, 4> dnaBases = {'A', 'C', 'G', 'T'};

  for (const auto& altBase : dnaBases) {
    if (altBase != seq[base_pos]) {
      auto altSeq = std::string(seq.cbegin(), seq.cend());
      altSeq[base_pos] = altBase;

      result.emplace_back(utils::RevComp(altSeq));
      result.emplace_back(std::move(altSeq));
    }
  }

  return result;
}
}  // namespace lancet2

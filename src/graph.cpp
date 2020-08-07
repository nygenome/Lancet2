#include "lancet/graph.h"

#include <algorithm>
#include <cmath>
#include <deque>

#include "absl/container/btree_set.h"
#include "lancet/align.h"
#include "lancet/assert_macro.h"
#include "lancet/canonical_kmers.h"
#include "lancet/dot_serializer.h"
#include "lancet/edmond_karp.h"
#include "lancet/logger.h"
#include "lancet/node_neighbour.h"
#include "lancet/tandem_repeat.h"
#include "lancet/timer.h"
#include "lancet/utils.h"
#include "lancet/variant.h"

namespace lancet {
Graph::Graph(std::shared_ptr<const RefWindow> w, Graph::NodeContainer&& data, double avg_cov, std::size_t k,
             std::shared_ptr<const CliParams> p)
    : window(std::move(w)), avgSampleCov(avg_cov), kmerSize(k), params(std::move(p)), nodesMap(std::move(data)) {}

void Graph::ProcessGraph(RefInfos&& ref_infos, const std::shared_ptr<VariantStore>& store) {
  Timer timer;
  const auto windowId = window->ToRegionString();
  DebugLog("Starting to process graph for %s with %d nodes", windowId, nodesMap.size());

  RemoveLowCovNodes(0);
  nodesMap.rehash(0);
  MarkConnectedComponents();

  for (const auto& comp : componentsInfo) {
    const auto markResult = MarkSourceSink(comp.ID);
    if (!markResult.foundSourceAndSink) continue;
    DebugLog("Marked source and sink in component%d (%d nodes) for %s", comp.ID, comp.numNodes, windowId);

    if (HasCycle()) {
      shouldIncrementK = true;
      DebugLog("Found graph cycle in component%d for %s with K=%d", comp.ID, windowId, kmerSize);
      return;
    }

    if (!params->outGraphsDir.empty()) WriteDot(comp.ID, "before_pruning");
    CompressGraph(comp.ID);
    RemoveLowCovNodes(comp.ID);
    CompressGraph(comp.ID);
    RemoveTips(comp.ID);
    RemoveShortLinks(comp.ID);
    nodesMap.rehash(0);
    if (!params->outGraphsDir.empty()) WriteDot(comp.ID, "after_pruning");

    if (HasCycle()) {
      shouldIncrementK = true;
      DebugLog("Found graph cycle in component%d for %s with K=%d", comp.ID, windowId, kmerSize);
      return;
    }

    LANCET_ASSERT(HasRefEnds(nodesMap, markResult, window->SeqView(), kmerSize));  // NOLINT

    std::size_t numPaths = 0;
    std::size_t numVariants = 0;
    const auto clampedRefInfos = ClampToSourceSink(ref_infos, markResult);
    const auto maxPathLength = RefAnchorLen(markResult) + static_cast<std::size_t>(params->maxIndelLength);

    EdmondKarpMaxFlow flow(&nodesMap, kmerSize, maxPathLength, params->graphTraversalLimit, params->tenxMode);
    std::vector<PathNodeIds> perPathTouches;
    auto pathPtr = flow.NextPath();

    while (pathPtr != nullptr) {
      numPaths++;
      if (!params->outGraphsDir.empty()) perPathTouches.emplace_back(pathPtr->TouchedEdgeIDs());

      if (utils::HasAlmostRepeatKmer(pathPtr->SeqView(), kmerSize, params->maxRptMismatch)) {
        DebugLog("Found repeat %d-mer in path%d of component%d for %s", kmerSize, numPaths, comp.ID, windowId);
        shouldIncrementK = true;
        return;
      }

      const auto transcripts = ProcessPath(*pathPtr, clampedRefInfos, markResult);
      numVariants += AddTranscripts(absl::MakeConstSpan(transcripts), store);
      pathPtr = flow.NextPath();
    }

    if (numPaths == 0) DebugLog("No path found in component%d for %s with K=%d", comp.ID, windowId, kmerSize);
    if (numVariants > 0) DebugLog("Added %d variant(s) from %s using K=%d", numVariants, windowId, kmerSize);
    if (!params->outGraphsDir.empty()) WriteDot(comp.ID, absl::MakeConstSpan(perPathTouches));
  }

  DebugLog("Done processing graph for %s | Runtime=%s", windowId, timer.HumanRuntime());
}

void Graph::MarkConnectedComponents() {
  std::size_t currentComponent = 0;
  componentsInfo.clear();

#ifndef NDEBUG
  static const auto isUnassigned = [](NodeContainer::const_reference p) { return p.second->ComponentID == 0; };
#endif
  LANCET_ASSERT(std::count_if(nodesMap.cbegin(), nodesMap.cend(), isUnassigned) == nodesMap.size());  // NOLINT

  for (NodeContainer::reference p : nodesMap) {
    if (p.second->ComponentID != 0) continue;

    currentComponent++;
    componentsInfo.emplace_back(ComponentInfo{currentComponent, 0});

    std::deque<Node*> connectedNodes;
    connectedNodes.push_back(p.second.get());

    while (!connectedNodes.empty()) {
      auto* currNode = connectedNodes.front();
      LANCET_ASSERT(currNode != nullptr);  // NOLINT

      if (currNode->ComponentID != 0) {
        connectedNodes.pop_front();
        continue;
      }

      currNode->ComponentID = currentComponent;
      componentsInfo[currentComponent - 1].numNodes += 1;
      for (const Edge& e : *currNode) {
        const auto neighbourItr = nodesMap.find(e.DestinationID());
        if (neighbourItr == nodesMap.end()) continue;
        LANCET_ASSERT(neighbourItr->second != nullptr);  // NOLINT
        connectedNodes.push_back(neighbourItr->second.get());
      }

      connectedNodes.pop_front();
    }
  }

  LANCET_ASSERT(std::count_if(nodesMap.cbegin(), nodesMap.cend(), isUnassigned) == 0);  // NOLINT
  DebugLog("Marked %d components in graph for %s", componentsInfo.size(), window->ToRegionString());
}

auto Graph::MarkSourceSink(std::size_t comp_id) -> Graph::MarkSourceSinkResult {
  auto refMerIDs = CanonicalKmerHashes(window->SeqView(), kmerSize);
  const auto srcResult = FindRefEnd(GraphEnd::SOURCE, comp_id, absl::MakeConstSpan(refMerIDs));
  if (!srcResult.foundEnd) return {false, 0, 0};

  std::reverse(refMerIDs.begin(), refMerIDs.end());
  const auto snkResult = FindRefEnd(GraphEnd::SINK, comp_id, absl::MakeConstSpan(refMerIDs));
  if (!snkResult.foundEnd) return {false, 0, 0};
  if (srcResult.nodeId == snkResult.nodeId) return {false, 0, 0};

  ResetSourceSink(nodesMap, comp_id);
  auto fauxSrcItr = nodesMap.find(MOCK_SOURCE_ID);
  auto fauxSnkItr = nodesMap.find(MOCK_SINK_ID);
  auto dataSrcItr = nodesMap.find(srcResult.nodeId);
  auto dataSnkItr = nodesMap.find(snkResult.nodeId);
  LANCET_ASSERT(fauxSrcItr != nodesMap.end());  // NOLINT
  LANCET_ASSERT(fauxSnkItr != nodesMap.end());  // NOLINT
  LANCET_ASSERT(dataSrcItr != nodesMap.end());  // NOLINT
  LANCET_ASSERT(dataSnkItr != nodesMap.end());  // NOLINT

  DisconnectEdges(dataSrcItr, nodesMap, ReverseStrand(dataSrcItr->second->Orientation()));
  DisconnectEdges(dataSnkItr, nodesMap, dataSnkItr->second->Orientation());

  const auto fauxSrcToDataSrcKind = MakeEdgeKind(Strand::FWD, dataSrcItr->second->Orientation());
  fauxSrcItr->second->EmplaceEdge(dataSrcItr->first, fauxSrcToDataSrcKind);
  dataSrcItr->second->EmplaceEdge(MOCK_SOURCE_ID, ReverseEdgeKind(fauxSrcToDataSrcKind));

  const auto fauxSnkToDataSnkKind = MakeEdgeKind(Strand::FWD, dataSnkItr->second->Orientation());
  fauxSnkItr->second->EmplaceEdge(dataSnkItr->first, fauxSnkToDataSnkKind);
  dataSnkItr->second->EmplaceEdge(MOCK_SINK_ID, ReverseEdgeKind(fauxSnkToDataSnkKind));

  const auto startBaseIdx = srcResult.refMerIdx;
  const auto endBaseIdx = snkResult.refMerIdx + dataSnkItr->second->Length();

  LANCET_ASSERT((window->SeqView().substr(startBaseIdx, kmerSize) == dataSrcItr->second->SeqView() ||
                 window->SeqView().substr(startBaseIdx, kmerSize) == dataSrcItr->second->FwdSeq()) &&
                (window->SeqView().substr(endBaseIdx - kmerSize, kmerSize) == dataSnkItr->second->SeqView() ||
                 window->SeqView().substr(endBaseIdx - kmerSize, kmerSize) == dataSnkItr->second->FwdSeq()));
  return MarkSourceSinkResult{true, startBaseIdx, endBaseIdx};
}

auto Graph::RemoveLowCovNodes(std::size_t comp_id) -> bool {
  // minNodeCov -> minimum coverage required for each node.
  // minWindowCov -> avg window coverage * MIN_NODE_COV_RATIO for each node
  const auto minWindowCov = static_cast<std::uint16_t>(std::ceil(params->minCovRatio * avgSampleCov));
  const auto minReqCov = std::max(static_cast<std::uint16_t>(params->minNodeCov), minWindowCov);

  std::vector<NodeIdentifier> nodesToRemove;
  std::for_each(nodesMap.cbegin(), nodesMap.cend(),
                [&nodesToRemove, &comp_id, &minReqCov](NodeContainer::const_reference p) {
                  if (p.second->IsMockNode() || p.second->ComponentID != comp_id) return;

                  const auto isNormalSingleton = p.second->SampleCount(SampleLabel::NORMAL) == 1;
                  const auto isTumorSingleton = p.second->SampleCount(SampleLabel::TUMOR) == 1;
                  if ((isNormalSingleton && isTumorSingleton) || p.second->MinSampleBaseCov() <= minReqCov) {
                    nodesToRemove.emplace_back(p.first);
                  }
                });

  if (!nodesToRemove.empty()) {
    DebugLog("Removing %d (%.2f%%) low cov nodes in component%d for %s", nodesToRemove.size(),
             100.0 * (static_cast<double>(nodesToRemove.size()) / static_cast<double>(nodesMap.size())), comp_id,
             window->ToRegionString());

    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
  }

  return !nodesToRemove.empty();
}

auto Graph::CompressGraph(std::size_t comp_id) -> bool {
  absl::flat_hash_set<NodeIdentifier> nodesToRemove;
  for (NodeContainer::const_reference p : nodesMap) {
    if (p.second->ComponentID != comp_id || p.second->IsMockNode()) continue;
    if (nodesToRemove.find(p.first) != nodesToRemove.end()) continue;
    CompressNode(p.first, FindCompressibleNeighbours(p.first), &nodesToRemove);
  }

  if (!nodesToRemove.empty()) {
    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
    DebugLog("Compressed %d nodes in component%d for %s", nodesToRemove.size(), comp_id, window->ToRegionString());
  }

  return !nodesToRemove.empty();
}

auto Graph::RemoveTips(std::size_t comp_id) -> bool {
  std::size_t totalTips = 0;
  std::size_t numTips = 0;
  const auto currK = kmerSize;
  const auto minTipLen = static_cast<std::size_t>(params->minGraphTipLength);
  const auto windowId = window->ToRegionString();

  // remove tips and compress at least once. compression after tip removal
  // can produce new tips in the graph, so continue until there are no tips
  do {
    std::vector<NodeIdentifier> nodesToRemove;
    std::for_each(nodesMap.cbegin(), nodesMap.cend(),
                  [&nodesToRemove, &comp_id, &currK, &minTipLen](NodeContainer::const_reference p) {
                    if (p.second->IsMockNode() || p.second->ComponentID == comp_id) return;
                    if (p.second->NumEdges() <= 1 && (p.second->Length() - currK + 1) < minTipLen) {
                      nodesToRemove.emplace_back(p.first);
                    }
                  });

    numTips = nodesToRemove.size();
    if (!nodesToRemove.empty()) {
      totalTips += numTips;
      RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
      CompressGraph(comp_id);
    }
  } while (numTips > 0);

  if (totalTips > 0) DebugLog("Removed %d tips in component%d for %s", totalTips, comp_id, windowId);
  return totalTips > 0;
}

auto Graph::RemoveShortLinks(std::size_t comp_id) -> bool {
  const auto currK = kmerSize;
  const auto minLinkLen = static_cast<std::size_t>(std::floor(static_cast<double>(kmerSize) / 2.0));
  const auto minReqCov = std::floor(std::sqrt(avgSampleCov));
  const TandemRepeatParams tandemParams{params->maxSTRUnitLength, params->minSTRUnits, params->minSTRLen,
                                        params->maxSTRDist};

  std::vector<NodeIdentifier> nodesToRemove;
  std::for_each(
      nodesMap.cbegin(), nodesMap.cend(),
      [&nodesToRemove, &comp_id, &currK, &minLinkLen, &minReqCov, &tandemParams](NodeContainer::const_reference p) {
        if (p.second->IsMockNode() || p.second->ComponentID == comp_id) return;

        const auto nodeDegree = p.second->NumEdges();
        const auto uniqSeqLen = p.second->Length() - currK + 1;
        const auto minRawCov = static_cast<float>(p.second->MinSampleBaseCov());
        if (nodeDegree >= 2 && uniqSeqLen < minLinkLen && minRawCov <= minReqCov) {
          const auto trQuery = FindTandemRepeat(p.second->SeqView(), currK - 1, tandemParams);
          // do not remove short-links within STRs: small bubbles are normal in STRs
          if (!trQuery.foundSTR) nodesToRemove.emplace_back(p.first);
        }
      });

  if (!nodesToRemove.empty()) {
    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
    DebugLog("Removed %d short links in component%d for %s", nodesToRemove.size(), comp_id, window->ToRegionString());
    CompressGraph(comp_id);
  }

  return !nodesToRemove.empty();
}

auto Graph::HasCycle() const -> bool {
  absl::flat_hash_set<NodeIdentifier> touchedIDs;
  return HasCycle(MOCK_SOURCE_ID, Strand::FWD, &touchedIDs) || HasCycle(MOCK_SOURCE_ID, Strand::REV, &touchedIDs);
}

auto Graph::ProcessPath(const Path& path, const RefInfos& ref_infos, const MarkSourceSinkResult& end_info) const
    -> Graph::TranscriptList {
  const auto pathSeq = path.SeqView();
  const auto refAnchorSeq = window->SeqView().substr(end_info.startOffset, RefAnchorLen(end_info));
  if (pathSeq == refAnchorSeq) return {};

  // check that reference seq length and reference data lengths are same
  LANCET_ASSERT(refAnchorSeq.length() == ref_infos[0].length());  // NOLINT
  LANCET_ASSERT(refAnchorSeq.length() == ref_infos[1].length());  // NOLINT

  AlignedSequences rawAlignedSeqs;  // need to create this because `goto`
  auto aligned = AlignedSequencesView{refAnchorSeq, pathSeq};
  if (utils::HammingDistWithin(refAnchorSeq, pathSeq, 5)) goto SkipLocalAlignment;  // NOLINT

  rawAlignedSeqs = Align(refAnchorSeq, pathSeq);
  aligned.ref = rawAlignedSeqs.ref;
  aligned.qry = rawAlignedSeqs.qry;

SkipLocalAlignment:
  LANCET_ASSERT(aligned.ref.length() == aligned.qry.length());  // NOLINT
  TrimEndGaps(&aligned);

  // 0-based reference anchor position in absolute chromosome coordinates
  const auto anchorStartGenome = static_cast<std::size_t>(window->StartPosition0()) + end_info.startOffset;
  std::size_t refIdx = 0;   // 0-based coordinate
  std::size_t refPos = 0;   // 1-based coordinate
  std::size_t pathPos = 0;  // 1-based coordinate

  auto code = TranscriptCode::REF_MATCH;
  TranscriptCode prevCode;
  TranscriptOffsets tmpOffsets;
  TranscriptBases tmpBases;

  std::vector<Transcript> transcripts;
  for (std::size_t idx = 0; idx < aligned.ref.length(); ++idx) {
    prevCode = code;

    if (aligned.ref[idx] == ALIGN_GAP) {
      code = TranscriptCode::INSERTION;
      refIdx = refPos;  // save variant position in reference before increment
      ++pathPos;
    } else if (aligned.qry[idx] == ALIGN_GAP) {
      code = TranscriptCode::DELETION;
      refIdx = refPos;  // save variant position in reference before increment
      ++refPos;
    } else {
      code = aligned.ref[idx] == aligned.qry[idx] ? TranscriptCode::REF_MATCH : TranscriptCode::SNV;
      refIdx = refPos;  // save variant position in reference before increment
      ++refPos;
      ++pathPos;
    }

    if (code == TranscriptCode::REF_MATCH) continue;

    const auto pathIdx = pathPos - 1;                          // 0-based index into the path sequence
    const auto genomeRefPos = anchorStartGenome + refIdx + 1;  // 1-based genome position

    const auto* spanner = path.FindSpanningNode(pathPos, kmerSize);
    LANCET_ASSERT(spanner != nullptr);  // NOLINT
    const auto withinTumorNode = spanner->LabelRatio(KmerLabel::TUMOR) >= 0.8;

    // compute previous base to the current event for both
    // ref and path sequence. [required for VCF output format]
    auto prevRefIdx = idx - 1;
    auto prevPathIdx = idx - 1;

    // must always be true because we force the ref-path alignment to always align at source and sink.
    LANCET_ASSERT(idx > 0);  // NOLINT
    while (aligned.ref[prevRefIdx] != 'A' && aligned.ref[prevRefIdx] != 'C' && aligned.ref[prevRefIdx] != 'G' &&
           aligned.ref[prevRefIdx] != 'T') {
      --prevRefIdx;
    }
    while (aligned.qry[prevPathIdx] != 'A' && aligned.qry[prevPathIdx] != 'C' && aligned.qry[prevPathIdx] != 'G' &&
           aligned.qry[prevPathIdx] != 'T') {
      --prevPathIdx;
    }

    LANCET_ASSERT(pathIdx < path.Length());         // NOLINT
    LANCET_ASSERT(refIdx < ref_infos[0].length());  // NOLINT
    LANCET_ASSERT(refIdx < ref_infos[1].length());  // NOLINT

    // create new transcript if we are sure that we can't extend a previous event
    if (transcripts.empty() || prevCode == TranscriptCode::REF_MATCH) {
      tmpOffsets.refStart = refIdx;
      tmpOffsets.altStart = pathIdx;
      tmpOffsets.refEnd = refIdx + 1;
      tmpOffsets.altEnd = pathIdx + 1;

      tmpBases.refBase = aligned.ref[idx];
      tmpBases.altBase = aligned.qry[idx];
      tmpBases.prevRefBase = aligned.ref[prevRefIdx];
      tmpBases.prevAltBase = aligned.qry[prevPathIdx];

      std::array<SampleCov, 2> sampleCovs{
          SampleCov(ref_infos[0].at(refIdx), path.HpCovAt(SampleLabel::NORMAL, pathIdx)),
          SampleCov(ref_infos[1].at(refIdx), path.HpCovAt(SampleLabel::TUMOR, pathIdx))};

      const auto chromName = window->Chromosome();
      transcripts.emplace_back(chromName, genomeRefPos, code, tmpOffsets, tmpBases, sampleCovs, withinTumorNode);
      continue;
    }

    // extend transcript from previous event
    Transcript& tr = transcripts[transcripts.size() - 1];
    const auto sameTranscriptCode = tr.code() == code;

    if (withinTumorNode && !tr.is_somatic()) tr.set_somatic_status(true);
    tr.add_ref_base(aligned.ref[idx]).add_alt_base(aligned.qry[idx]);
    if (code == TranscriptCode::INSERTION || code == TranscriptCode::SNV) tr.set_alt_end_offset(pathIdx + 1);
    if (code == TranscriptCode::DELETION || code == TranscriptCode::SNV) tr.set_ref_end_offset(refIdx + 1);

    // extend existing insertion, if possible
    if (sameTranscriptCode && code == TranscriptCode::INSERTION && tr.position() == genomeRefPos) {
      tr.add_coverage(SampleLabel::TUMOR, Allele::ALT, path.HpCovAt(SampleLabel::TUMOR, pathIdx))
          .add_coverage(SampleLabel::NORMAL, Allele::ALT, path.HpCovAt(SampleLabel::NORMAL, pathIdx));
      continue;
    }

    // extend existing deletion, if possible
    const auto deletedRefLen = tr.alt_seq().length();
    if (sameTranscriptCode && code == TranscriptCode::DELETION && (tr.position() + deletedRefLen) == genomeRefPos) {
      tr.add_coverage(SampleLabel::NORMAL, Allele::REF, ref_infos[0].at(refIdx))
          .add_coverage(SampleLabel::TUMOR, Allele::REF, ref_infos[1].at(refIdx));
      continue;
    }

    // extend into MNP or complex event
    // If current code is SNV & previous code is SNV, extend into MNP (also complex event for now)
    // If current code is SNV & previous code is not SNV, extend into complex event
    tr.set_code(TranscriptCode::COMPLEX)
        .add_coverage(SampleLabel::NORMAL, Allele::REF, ref_infos[0].at(refIdx))
        .add_coverage(SampleLabel::TUMOR, Allele::REF, ref_infos[1].at(refIdx))
        .add_coverage(SampleLabel::TUMOR, Allele::ALT, path.HpCovAt(SampleLabel::TUMOR, pathIdx))
        .add_coverage(SampleLabel::NORMAL, Allele::ALT, path.HpCovAt(SampleLabel::NORMAL, pathIdx));
  }

  // If alignment left shifts the InDel, reference and path coverages can get out of sync.
  // Add coverage for k-1 bases after reference and path ends to fix this.
  const auto k = kmerSize;  // tmp for lambda
  const TandemRepeatParams tandemParams{params->maxSTRUnitLength, params->minSTRUnits, params->minSTRLen,
                                        params->maxSTRDist};

  std::for_each(
      transcripts.begin(), transcripts.end(), [&path, &ref_infos, &k, &tandemParams, &pathSeq](Transcript& transcript) {
        transcript.add_str_result(FindTandemRepeat(pathSeq, transcript.alt_start_offset(), tandemParams));

        if (transcript.code() == TranscriptCode::REF_MATCH || transcript.code() == TranscriptCode::SNV) return;

        for (std::size_t pos = 0; pos <= k; pos++) {
          const auto currPathIdx = transcript.alt_end_offset() + pos;
          const auto currRefIdx = transcript.ref_end_offset() + pos;

          const auto* spanner = path.FindSpanningNode(currPathIdx, k);
          LANCET_ASSERT(spanner != nullptr);  // NOLINT
          constexpr double minRatioForSomatic = 0.8;
          if (spanner->LabelRatio(KmerLabel::TUMOR) >= minRatioForSomatic) transcript.set_somatic_status(true);

          if (currRefIdx < ref_infos[0].length() && currRefIdx < ref_infos[1].length()) {
            transcript.add_coverage(SampleLabel::NORMAL, Allele::REF, ref_infos[0].at(currRefIdx))
                .add_coverage(SampleLabel::TUMOR, Allele::REF, ref_infos[1].at(currRefIdx));
          }

          if (currPathIdx >= path.Length()) continue;
          transcript.add_coverage(SampleLabel::TUMOR, Allele::ALT, path.HpCovAt(SampleLabel::TUMOR, currPathIdx))
              .add_coverage(SampleLabel::NORMAL, Allele::ALT, path.HpCovAt(SampleLabel::NORMAL, currPathIdx));
        }
      });

  return transcripts;
}

auto Graph::AddTranscripts(absl::Span<const Transcript> transcripts, const std::shared_ptr<VariantStore>& store) const
    -> std::size_t {
  std::size_t numVariantsAdded = 0;
  const auto windowIdx = window->WindowIndex();

  for (const auto& transcript : transcripts) {
    if (!transcript.has_alt_coverage()) continue;
    Variant variant(transcript, kmerSize);
    if (variant.ComputeState() == VariantState::NONE) continue;
    if (store->AddVariant(windowIdx, std::move(variant))) numVariantsAdded++;
  }

  return numVariantsAdded;
}

void Graph::WriteDot(std::size_t comp_id, const std::string& suffix) const {
  const DotSerializer ds(this);
  ds.write_component(comp_id, suffix);
}

void Graph::WriteDot(std::size_t comp_id, absl::Span<const PathNodeIds> flow_paths) const {
  const DotSerializer ds(this);
  ds.write_component(comp_id, flow_paths);
}

void Graph::EraseNode(NodeIterator itr) {
  if (itr == nodesMap.end() || itr->second->IsMockNode()) return;

  // remove edges associated with the to be removed nodes first, then remove the node
  for (const Edge& e : *itr->second) {
    auto neighbourItr = nodesMap.find(e.DestinationID());
    if (neighbourItr == nodesMap.end()) continue;
    neighbourItr->second->EraseEdge(itr->first, ReverseEdgeKind(e.Kind()));
  }

  nodesMap.erase(itr);
}

void Graph::EraseNode(NodeIdentifier node_id) { return EraseNode(nodesMap.find(node_id)); }

auto Graph::FindRefEnd(GraphEnd k, std::size_t comp_id, absl::Span<const NodeIdentifier> ref_mer_hashes) const
    -> Graph::RefEndResult {
  // find node in component for the first/last occurring ref kmer id.
  const auto minEndCov = static_cast<std::uint16_t>(params->minAnchorCov);
  const auto numRefMers = static_cast<std::int64_t>(ref_mer_hashes.size());

  for (auto refIdx = 0; refIdx < numRefMers; ++refIdx) {
    const auto merIndex = static_cast<std::size_t>(refIdx);
    const auto itr = nodesMap.find(ref_mer_hashes[merIndex]);
    if (itr == nodesMap.end()) continue;

    LANCET_ASSERT(itr->second != nullptr && !itr->second->IsMockNode());  // NOLINT
    if (itr->second->ComponentID != comp_id || itr->second->TotalSampleCount() < minEndCov) continue;

    const auto resultMerIdx = k == GraphEnd::SOURCE ? merIndex : (numRefMers - merIndex - 1);
    return RefEndResult{itr->first, resultMerIdx, true};
  }

  return {0, 0, false};
}

auto Graph::FindCompressibleNeighbours(NodeIdentifier src_id) const -> absl::btree_set<NodeNeighbour> {
  if (src_id == MOCK_SOURCE_ID || src_id == MOCK_SINK_ID) return {};

  const auto srcItr = nodesMap.find(src_id);
  LANCET_ASSERT(srcItr != nodesMap.end() && srcItr->second != nullptr);  // NOLINT
  const auto srcNeighbours = srcItr->second->FindMergeableNeighbours();
  if (srcNeighbours.empty()) return {};

  absl::btree_set<NodeNeighbour> results;
  for (const auto& srcNbour : srcNeighbours) {
    const auto buddyItr = nodesMap.find(srcNbour.buddyId);
    LANCET_ASSERT(buddyItr != nodesMap.end() && buddyItr->second != nullptr);  // NOLINT
    const auto buddysNeighbours = buddyItr->second->FindMergeableNeighbours();
    const auto areMututalBuddies = std::any_of(buddysNeighbours.cbegin(), buddysNeighbours.cend(),
                                               [&src_id](const NodeNeighbour& n) { return n.buddyId == src_id; });

    const auto mergeDir = SourceStrand(srcNbour.edgeKind) == Strand::FWD ? BuddyPosition::FRONT : BuddyPosition::BACK;
    const auto canMergeWithBuddySeq = srcItr->second->CanMerge(*buddyItr->second, mergeDir, kmerSize);
    if (areMututalBuddies && canMergeWithBuddySeq) results.emplace(srcNbour);
  }

  return results;
}

void Graph::CompressNode(NodeIdentifier src_id, const absl::btree_set<NodeNeighbour>& buddies,
                         absl::flat_hash_set<NodeIdentifier>* compressed) const {
  if (buddies.empty() || buddies.size() > 2) return;

  const auto srcItr = nodesMap.find(src_id);
  LANCET_ASSERT(srcItr != nodesMap.end());  // NOLINT

  absl::btree_set<NodeNeighbour> remaining;
  std::for_each(buddies.cbegin(), buddies.cend(), [&compressed, &remaining](const NodeNeighbour& n) {
    if (compressed->find(n.buddyId) == compressed->end()) remaining.emplace(n);
  });

  while (!remaining.empty() && remaining.size() <= 2) {
    const auto srcToBuddy = *remaining.cbegin();
    LANCET_ASSERT(compressed->find(srcToBuddy.buddyId) == compressed->end());  // NOLINT
    const auto buddyItr = nodesMap.find(srcToBuddy.buddyId);
    LANCET_ASSERT(buddyItr != nodesMap.end());  // NOLINT

    const auto mergeDir = SourceStrand(srcToBuddy.edgeKind) == Strand::FWD ? BuddyPosition::FRONT : BuddyPosition::BACK;
    if (!srcItr->second->CanMerge(*buddyItr->second, mergeDir, kmerSize)) {
      remaining.erase(remaining.begin());
      continue;
    }

    srcItr->second->MergeBuddy(*buddyItr->second, mergeDir, kmerSize);
    srcItr->second->EraseEdge(srcToBuddy.buddyId);
    compressed->emplace(srcToBuddy.buddyId);

    const auto srcBuddyDiffStrands = SourceStrand(srcToBuddy.edgeKind) != DestStrand(srcToBuddy.edgeKind);
    for (const Edge& buddyE : *buddyItr->second) {
      const auto buddyNeighbourId = buddyE.DestinationID();
      if (buddyNeighbourId == src_id) continue;

      auto buddysNeighbourItr = nodesMap.find(buddyNeighbourId);
      LANCET_ASSERT(buddysNeighbourItr != nodesMap.end());  // NOLINT

      const auto srcLinkStrand = srcBuddyDiffStrands ? ReverseStrand(buddyE.SrcDirection()) : buddyE.SrcDirection();
      const auto resultKind = MakeEdgeKind(srcLinkStrand, buddyE.DstDirection());

      if (buddyNeighbourId == srcToBuddy.buddyId) {
        srcItr->second->EmplaceEdge(srcItr->first, resultKind);
        continue;
      }

      srcItr->second->EmplaceEdge(buddyNeighbourId, resultKind);
      buddysNeighbourItr->second->EraseEdge(srcToBuddy.buddyId);
      buddysNeighbourItr->second->EmplaceEdge(src_id, ReverseEdgeKind(resultKind));
    }

    remaining.erase(remaining.begin());
    const auto newNeighbours = FindCompressibleNeighbours(srcItr->first);
    std::for_each(newNeighbours.cbegin(), newNeighbours.cend(), [&compressed, &remaining](const NodeNeighbour& n) {
      if (compressed->find(n.buddyId) == compressed->end()) remaining.emplace(n);
    });
  }
}

auto Graph::HasCycle(NodeIdentifier node_id, Strand direction, absl::flat_hash_set<NodeIdentifier>* touched) const
    -> bool {
  const auto itr = nodesMap.find(node_id);
  if (itr == nodesMap.end()) return false;

  touched->insert(node_id);
  for (const Edge& e : *itr->second) {
    const auto neighbourId = e.DestinationID();
    if (neighbourId == MOCK_SOURCE_ID || neighbourId == MOCK_SINK_ID || e.SrcDirection() != direction) continue;
    if (touched->find(neighbourId) == touched->end()) return HasCycle(neighbourId, e.DstDirection(), touched);
    touched->erase(node_id);
    return true;
  }

  touched->erase(node_id);
  return false;
}

auto Graph::ClampToSourceSink(const RefInfos& refs, const MarkSourceSinkResult& ends) -> Graph::RefInfos {
  const auto length = ends.endOffset - ends.startOffset;
  return std::array<absl::Span<const BaseHpCov>, 2>{refs[0].subspan(ends.startOffset, length),
                                                    refs[1].subspan(ends.startOffset, length)};
}

void Graph::ResetSourceSink(const NodeContainer& nc, std::size_t current_component) {
  const auto fauxSrcItr = nc.find(MOCK_SOURCE_ID);
  const auto fauxSnkItr = nc.find(MOCK_SINK_ID);
  LANCET_ASSERT(fauxSrcItr != nc.end());  // NOLINT
  LANCET_ASSERT(fauxSnkItr != nc.end());  // NOLINT

  fauxSrcItr->second->ComponentID = current_component;
  fauxSnkItr->second->ComponentID = current_component;
  fauxSrcItr->second->ClearEdges();
  fauxSnkItr->second->ClearEdges();
}

void Graph::DisconnectEdges(NodeIterator itr, const NodeContainer& nc, Strand direction) {
  LANCET_ASSERT(itr != nc.end());  // NOLINT

  for (const Edge& e : *itr->second) {
    if (e.DestinationID() == itr->first || e.SrcDirection() != direction) continue;
    auto neighbourItr = nc.find(e.DestinationID());
    if (neighbourItr == nc.end()) continue;

    itr->second->EraseEdge(neighbourItr->first, e.Kind());
    neighbourItr->second->EraseEdge(itr->first);
  }
}

auto Graph::HasRefEnds(const NodeContainer& nc, const MarkSourceSinkResult& ends, std::string_view ref, std::size_t k)
    -> bool {
  const auto fauxSrcItr = nc.find(MOCK_SOURCE_ID);
  const auto fauxSnkItr = nc.find(MOCK_SINK_ID);

  LANCET_ASSERT(fauxSrcItr != nc.end() && fauxSnkItr != nc.end());
  if (fauxSrcItr == nc.end() || fauxSnkItr == nc.end()) return false;

  LANCET_ASSERT(fauxSrcItr->second->NumEdges() == 1 && fauxSnkItr->second->NumEdges() == 1);
  const auto dataSrcItr = nc.find((*fauxSrcItr->second->begin()).DestinationID());
  const auto dataSnkItr = nc.find((*fauxSnkItr->second->begin()).DestinationID());

  LANCET_ASSERT(dataSrcItr != nc.end() && dataSnkItr != nc.end());
  if (dataSrcItr == nc.end() || dataSnkItr == nc.end()) return false;

  const auto snkLen = dataSnkItr->second->Length();
  return (ref.substr(ends.startOffset, k) == dataSrcItr->second->SeqView().substr(0, k) ||
          ref.substr(ends.startOffset, k) == dataSrcItr->second->FwdSeq().substr(0, k)) &&
         (ref.substr(ends.endOffset - k, k) == dataSnkItr->second->SeqView().substr(snkLen - k, k) ||
          ref.substr(ends.endOffset - k, k) == dataSnkItr->second->FwdSeq().substr(snkLen - k, k));
}
}  // namespace lancet

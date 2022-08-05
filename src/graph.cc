#include "lancet2/graph.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "absl/container/btree_set.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "lancet2/align.h"
#include "lancet2/assert_macro.h"
#include "lancet2/canonical_kmers.h"
#include "lancet2/dot_serializer.h"
#include "lancet2/edmond_karp.h"
#include "lancet2/kmer.h"
#include "lancet2/log_macros.h"
#include "lancet2/node_neighbour.h"
#include "lancet2/tandem_repeat.h"
#include "lancet2/timer.h"
#include "lancet2/utils.h"
#include "lancet2/variant.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
Graph::Graph(std::shared_ptr<const RefWindow> w, Graph::NodeContainer&& data, absl::Span<const ReadInfo> reads,
             double avg_cov, usize k, std::shared_ptr<const CliParams> p)
    : kmerSize(k), avgSampleCov(avg_cov), params(std::move(p)), window(std::move(w)), sampleReads(reads),
      nodesMap(std::move(data)) {}

void Graph::ProcessGraph(std::vector<Variant>* results) {
  Timer timer;
  const auto windowId = window->ToRegionString();
  LOG_DEBUG("Starting to process graph for {} with {} nodes", windowId, nodesMap.size());

  RemoveLowCovNodes(0);
  nodesMap.rehash(0);
  const auto componentsInfo = MarkConnectedComponents();

  for (const auto& comp : componentsInfo) {
    const auto markResult = MarkSourceSink(comp.ID);
    if (!markResult.foundSrcAndSnk) continue;
    LOG_DEBUG("Marked source and sink in component{} ({} nodes) for {}", comp.ID, comp.numNodes, windowId);

    if (HasCycle()) {
      shouldIncrementK = true;
      LOG_DEBUG("Found graph cycle in component{} for {} with K={}", comp.ID, windowId, kmerSize);
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
      LOG_DEBUG("Found graph cycle in component{} for {} with K={}", comp.ID, windowId, kmerSize);
      return;
    }

    usize numPaths = 0;
    const auto maxPathLength = RefAnchorLen(markResult) + static_cast<usize>(params->maxIndelLength);

    EdmondKarpMaxFlow flow(&nodesMap, kmerSize, maxPathLength, params->graphTraversalLimit, params->tenxMode);
    std::vector<PathNodeIds> perPathTouches;
    auto pathPtr = flow.NextPath();

    while (pathPtr != nullptr) {
      numPaths++;
      if (!params->outGraphsDir.empty()) perPathTouches.emplace_back(pathPtr->TouchedEdgeIDs());

      if (utils::HasAlmostRepeatKmer(pathPtr->GetSeqView(), kmerSize, params->maxRptMismatch)) {
        LOG_DEBUG("Found repeat {}-mer in path{} of component{} for {}", kmerSize, numPaths, comp.ID, windowId);
        shouldIncrementK = true;
        return;
      }

      ProcessPath(*pathPtr, markResult, results);
      if (!params->outGraphsDir.empty()) {
        WritePathFasta(pathPtr->GetSeqView(), comp.ID, numPaths);
      }

      pathPtr = flow.NextPath();
    }

    if (numPaths == 0) LOG_DEBUG("No path found in component{} for {} with K={}", comp.ID, windowId, kmerSize);
    if (!params->outGraphsDir.empty() && !perPathTouches.empty()) {
      WriteDot(comp.ID, absl::MakeConstSpan(perPathTouches));
    }
  }

  LOG_DEBUG("Done processing graph for {} | Runtime={}", windowId, timer.HumanRuntime());
}

auto Graph::MarkConnectedComponents() -> std::vector<ComponentInfo> {
  usize currentComponent = 0;
  std::vector<ComponentInfo> componentsInfo;

#ifndef NDEBUG
  static const auto isUnassigned = [](NodeContainer::const_reference p) { return p.second->ComponentID == 0; };
#endif
  LANCET_ASSERT(static_cast<usize>(std::count_if(nodesMap.cbegin(), nodesMap.cend(), isUnassigned)) ==
                nodesMap.size());  // NOLINT

  for (NodeContainer::reference p : nodesMap) {
    if (p.second->ComponentID != 0) continue;

    currentComponent++;
    componentsInfo.emplace_back(ComponentInfo{currentComponent, 0});

    std::deque<Node*> connectedNodes{};
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
        const auto neighbourItr = nodesMap.find(e.GetDstID());
        if (neighbourItr == nodesMap.end()) continue;
        LANCET_ASSERT(neighbourItr->second != nullptr);  // NOLINT
        connectedNodes.push_back(neighbourItr->second.get());
      }

      connectedNodes.pop_front();
    }
  }

  LANCET_ASSERT(std::count_if(nodesMap.cbegin(), nodesMap.cend(), isUnassigned) == 0);  // NOLINT
  LOG_DEBUG("Marked {} components in graph for {}", componentsInfo.size(), window->ToRegionString());

  return componentsInfo;
}

auto Graph::MarkSourceSink(usize comp_id) -> Graph::SrcSnkResult {
  const auto refseq = window->GetSeqView();
  auto refMerIDs = CanonicalKmerHashes(refseq, kmerSize);
  const auto srcResult = FindRefEnd(GraphEnd::SOURCE, comp_id, absl::MakeConstSpan(refMerIDs));
  if (!srcResult.foundEnd) return {false, 0, 0};

  std::reverse(refMerIDs.begin(), refMerIDs.end());
  const auto snkResult = FindRefEnd(GraphEnd::SINK, comp_id, absl::MakeConstSpan(refMerIDs));
  if (!snkResult.foundEnd || srcResult.nodeId == snkResult.nodeId) return {false, 0, 0};

  auto fauxSrcItr = nodesMap.find(MOCK_SOURCE_ID);
  LANCET_ASSERT(fauxSrcItr != nodesMap.end());  // NOLINT
  fauxSrcItr->second->ComponentID = comp_id;
  DisconnectEdgesTo(fauxSrcItr, nodesMap);
  fauxSrcItr->second->ClearEdges();

  auto fauxSnkItr = nodesMap.find(MOCK_SINK_ID);
  LANCET_ASSERT(fauxSnkItr != nodesMap.end());  // NOLINT
  fauxSnkItr->second->ComponentID = comp_id;
  DisconnectEdgesTo(fauxSnkItr, nodesMap);
  fauxSnkItr->second->ClearEdges();

  auto dataSrcItr = nodesMap.find(srcResult.nodeId);
  auto dataSnkItr = nodesMap.find(snkResult.nodeId);
  LANCET_ASSERT(dataSrcItr != nodesMap.end());  // NOLINT
  LANCET_ASSERT(dataSnkItr != nodesMap.end());  // NOLINT

  const auto fauxSrcToDataSrcKind = MakeEdgeKind(Strand::FWD, dataSrcItr->second->GetOrientation());
  fauxSrcItr->second->EmplaceEdge(dataSrcItr->first, fauxSrcToDataSrcKind);
  dataSrcItr->second->EmplaceEdge(MOCK_SOURCE_ID, ReverseEdgeKind(fauxSrcToDataSrcKind));

  const auto isDataSnkRev = dataSnkItr->second->GetOrientation() == Strand::REV;
  const auto fauxSnkToDataSnkKind = isDataSnkRev ? EdgeKind::FF : EdgeKind::RR;
  fauxSnkItr->second->EmplaceEdge(dataSnkItr->first, fauxSnkToDataSnkKind);
  dataSnkItr->second->EmplaceEdge(MOCK_SINK_ID, ReverseEdgeKind(fauxSnkToDataSnkKind));

  const auto startBaseIdx = srcResult.refMerIdx;
  const auto endBaseIdx = snkResult.refMerIdx + dataSnkItr->second->GetLength();

  LANCET_ASSERT(fauxSrcItr->second->NumEdges() == 1);
  LANCET_ASSERT(fauxSnkItr->second->NumEdges() == 1);
  LANCET_ASSERT((refseq.substr(startBaseIdx, kmerSize) == dataSrcItr->second->GetSeqView() ||
                 utils::RevComp(refseq.substr(startBaseIdx, kmerSize)) == dataSrcItr->second->GetSeqView()) &&
                (refseq.substr(endBaseIdx - kmerSize, kmerSize) == dataSnkItr->second->GetSeqView() ||
                 utils::RevComp(refseq.substr(endBaseIdx - kmerSize, kmerSize)) == dataSnkItr->second->GetSeqView()));

  return SrcSnkResult{true, startBaseIdx, endBaseIdx};
}

auto Graph::RemoveLowCovNodes(usize comp_id) -> bool {
  // minNodeCov -> minimum coverage required for each node.
  // minWindowCov -> avg window coverage * MIN_NODE_COV_RATIO for each node
  const auto minWindowCov = static_cast<u16>(std::ceil(params->minCovRatio * avgSampleCov));
  const auto minReqCov = std::max(static_cast<u16>(params->minNodeCov), minWindowCov);

  std::vector<NodeIdentifier> nodesToRemove{};
  std::for_each(nodesMap.cbegin(), nodesMap.cend(),
                [&nodesToRemove, &comp_id, &minReqCov](NodeContainer::const_reference p) {
                  if (p.second->IsMockNode() || p.second->ComponentID != comp_id) return;

                  const auto isNormalSingleton = p.second->SampleCount(SampleLabel::NORMAL) == 1;
                  const auto isTumorSingleton = p.second->SampleCount(SampleLabel::TUMOR) == 1;
                  const auto totalSampleCov = p.second->TotalSampleCount();
                  if ((isNormalSingleton && isTumorSingleton) || totalSampleCov <= minReqCov) {
                    nodesToRemove.emplace_back(p.first);
                  }
                });

  if (!nodesToRemove.empty()) {
    LOG_DEBUG("Removing {} ({:.2f}%) low cov nodes in component{} for {}", nodesToRemove.size(),
              100.0 * (static_cast<double>(nodesToRemove.size()) / static_cast<double>(nodesMap.size())), comp_id,
              window->ToRegionString());

    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
  }

  return !nodesToRemove.empty();
}

auto Graph::CompressGraph(usize comp_id) -> bool {
  absl::flat_hash_set<NodeIdentifier> nodesToRemove;
  for (NodeContainer::const_reference p : nodesMap) {
    if (p.second->ComponentID != comp_id || p.second->IsMockNode()) continue;
    if (nodesToRemove.find(p.first) != nodesToRemove.end()) continue;
    CompressNode(p.first, FindCompressibleNeighbours(p.first), &nodesToRemove);
  }

  if (!nodesToRemove.empty()) {
    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
    LOG_DEBUG("Compressed {} nodes in component{} for {}", nodesToRemove.size(), comp_id, window->ToRegionString());
  }

  return !nodesToRemove.empty();
}

auto Graph::RemoveTips(usize comp_id) -> bool {
  usize totalTips = 0;
  usize numTips = 0;
  const auto currK = kmerSize;
  const auto minTipLen = static_cast<usize>(params->minGraphTipLength);
  const auto windowId = window->ToRegionString();

  // remove tips and compress at least once. compression after tip removal
  // can produce new tips in the graph, so continue until there are no tips
  do {
    std::vector<NodeIdentifier> nodesToRemove;
    std::for_each(nodesMap.cbegin(), nodesMap.cend(),
                  [&nodesToRemove, &comp_id, &currK, &minTipLen](NodeContainer::const_reference p) {
                    if (p.second->IsMockNode() || p.second->ComponentID != comp_id) return;
                    if (p.second->NumEdges() <= 1 && (p.second->GetLength() - currK + 1) < minTipLen) {
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

  if (totalTips > 0) LOG_DEBUG("Removed {} tips in component{} for {}", totalTips, comp_id, windowId);
  return totalTips > 0;
}

auto Graph::RemoveShortLinks(usize comp_id) -> bool {
  const auto currK = kmerSize;
  const auto minLinkLen = static_cast<usize>(std::floor(static_cast<double>(kmerSize) / 2.0));
  const auto minReqCov = std::floor(std::sqrt(avgSampleCov));
  const TandemRepeatParams tandemParams{params->maxSTRUnitLength, params->minSTRUnits, params->minSTRLen,
                                        params->maxSTRDist};

  std::vector<NodeIdentifier> nodesToRemove;
  std::for_each(
      nodesMap.cbegin(), nodesMap.cend(),
      [&nodesToRemove, &comp_id, &currK, &minLinkLen, &minReqCov, &tandemParams](NodeContainer::const_reference p) {
        if (p.second->IsMockNode() || p.second->ComponentID == comp_id) return;

        const auto nodeDegree = p.second->NumEdges();
        const auto uniqSeqLen = p.second->GetLength() - currK + 1;
        const auto minRawCov = static_cast<float>(p.second->TotalSampleCount());
        if (nodeDegree >= 2 && uniqSeqLen < minLinkLen && minRawCov <= minReqCov) {
          const auto trQuery = FindTandemRepeat(p.second->GetSeqView(), currK - 1, tandemParams);
          // do not remove short-links within STRs: small bubbles are normal in STRs
          if (!trQuery.foundSTR) nodesToRemove.emplace_back(p.first);
        }
      });

  if (!nodesToRemove.empty()) {
    RemoveNodes(nodesToRemove.cbegin(), nodesToRemove.cend());
    LOG_DEBUG("Removed {} short links in component{} for {}", nodesToRemove.size(), comp_id, window->ToRegionString());
    CompressGraph(comp_id);
  }

  return !nodesToRemove.empty();
}

auto Graph::HasCycle() const -> bool {
  absl::flat_hash_set<NodeIdentifier> touchedIDs;
  return HasCycle(MOCK_SOURCE_ID, Strand::FWD, &touchedIDs) || HasCycle(MOCK_SOURCE_ID, Strand::REV, &touchedIDs);
}

void Graph::ProcessPath(const Path& path, const SrcSnkResult& einfo, std::vector<Variant>* variants) const {
  const auto pathSeq = path.GetSeqView();
  const auto refAnchorSeq = window->GetSeqView().substr(einfo.startOffset, RefAnchorLen(einfo));
  if (pathSeq == refAnchorSeq) return;

  AlnSeqs rawAlignedSeqs;  // need to create this because `goto`
  auto aligned = AlnSeqsView{refAnchorSeq, pathSeq};
  if (utils::HammingDistWithin(refAnchorSeq, pathSeq, 5)) goto SkipLocalAlignment;  // NOLINT

  try {
    rawAlignedSeqs = Align(refAnchorSeq, pathSeq);
  } catch (const std::exception& e) {
    LOG_TRACE("Error processing window {}: error aligning ref: {}, qry: {} | exception – {}", window->ToRegionString(),
              refAnchorSeq, pathSeq, e.what());
    return;
  }

  aligned.ref = rawAlignedSeqs.ref;
  aligned.qry = rawAlignedSeqs.qry;

SkipLocalAlignment:
  LANCET_ASSERT(aligned.ref.length() == aligned.qry.length());  // NOLINT
  const auto trimAl = TrimEndGaps(&aligned);

  // 0-based reference anchor position in absolute chromosome coordinates
  const auto anchorGenomeStart = static_cast<usize>(window->GetStartPos0()) + einfo.startOffset + trimAl.refStartTrim;
  usize refIdx = 0;   // 0-based coordinate
  usize pathIdx = 0;  // 0-based coordinate
  usize refPos = 0;   // 1-based coordinate
  usize pathPos = 0;  // 1-based coordinate

  auto code = TranscriptCode::REF_MATCH;
  TranscriptCode prevCode = TranscriptCode::REF_MATCH;
  Transcript::Offsets tmpOffsets;
  Transcript::Bases tmpBases;

  std::vector<Transcript> transcripts;
  for (usize idx = 0; idx < aligned.ref.length(); ++idx) {
    prevCode = code;

    if (aligned.ref[idx] == ALIGN_GAP) {
      code = TranscriptCode::INSERTION;
      refIdx = refPos;  // save variant position in reference before increment
      pathIdx = pathPos;
      ++pathPos;
    } else if (aligned.qry[idx] == ALIGN_GAP) {
      code = TranscriptCode::DELETION;
      refIdx = refPos;  // save variant position in reference before increment
      pathIdx = pathPos;
      ++refPos;
    } else {
      code = aligned.ref[idx] == aligned.qry[idx] ? TranscriptCode::REF_MATCH : TranscriptCode::SNV;
      refIdx = refPos;  // save variant position in reference before increment
      pathIdx = pathPos;
      ++refPos;
      ++pathPos;
    }

    if (code == TranscriptCode::REF_MATCH) continue;

    const auto genomeRefPos = anchorGenomeStart + refIdx + 1;  // 1-based genome position

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

    LANCET_ASSERT(pathIdx < path.GetLength());  // NOLINT

    // create new transcript if we are sure that we can't extend a previous event
    if (transcripts.empty() || prevCode == TranscriptCode::REF_MATCH) {
      tmpOffsets.refStart = refIdx + trimAl.refStartTrim;
      tmpOffsets.altStart = pathIdx + trimAl.qryStartTrim;
      tmpOffsets.refEnd = refIdx + 1 + trimAl.refStartTrim;
      tmpOffsets.altEnd = pathIdx + 1 + trimAl.qryStartTrim;

      tmpBases.refBase = aligned.ref[idx];
      tmpBases.altBase = aligned.qry[idx];
      tmpBases.prevRefBase = aligned.ref[prevRefIdx];
      tmpBases.prevAltBase = aligned.qry[prevPathIdx];

      const auto chromName = window->GetChromName();
      transcripts.emplace_back(chromName, genomeRefPos, code, tmpOffsets, tmpBases);
      continue;
    }

    // extend transcript from previous event
    Transcript& tr = transcripts.back();
    const auto sameTranscriptCode = tr.Code() == code;

    tr.AddRefBase(aligned.ref[idx]).AddAltBase(aligned.qry[idx]);
    if (code == TranscriptCode::INSERTION || code == TranscriptCode::SNV) tr.SetAltEndOffset(pathIdx + 1);
    if (code == TranscriptCode::DELETION || code == TranscriptCode::SNV) tr.SetRefEndOffset(refIdx + 1);

    // extend existing insertion, if possible
    if (sameTranscriptCode && code == TranscriptCode::INSERTION && tr.Position() == genomeRefPos) {
      continue;
    }

    // extend existing deletion, if possible
    const auto deletedRefLen = tr.RefSeq().length();
    if (sameTranscriptCode && code == TranscriptCode::DELETION && (tr.Position() + deletedRefLen - 1) == genomeRefPos) {
      continue;
    }

    // extend into MNP or complex event
    // If current code is SNV & previous code is SNV, extend into MNP (also complex event for now)
    // If current code is SNV & previous code is not SNV, extend into complex event
    tr.SetCode(TranscriptCode::COMPLEX);
  }

  const TandemRepeatParams trP{params->maxSTRUnitLength, params->minSTRUnits, params->minSTRLen, params->maxSTRDist};
  for (Transcript& T : transcripts) {
    T.AddSTRResult(FindTandemRepeat(pathSeq, T.AltStartOffset(), trP));

    // Populate ref and alt haplotype spanning the ref and alt alleles with
    // length atleast kmer size or ref/alt allele length whichever is larger
    T.Finalize();
    T.BuildHaplotypes(refAnchorSeq, pathSeq, kmerSize);
  }

  BuildVariants(absl::MakeConstSpan(transcripts), variants);
}

void Graph::WritePathFasta(std::string_view path_seq, usize comp_id, usize path_num) const {
  const auto outDir =
      params->outGraphsDir.empty() ? std::filesystem::current_path() : std::filesystem::path(params->outGraphsDir);
  const auto windowID = window->ToRegionString();

  const auto pathID = absl::StrFormat("%s_c%d_p%d", windowID, comp_id, path_num);
  const auto fName = absl::StrFormat("%s/%s.fasta", outDir, windowID);

  std::ofstream outStream(fName, std::ios::out | std::ios::app);
  outStream << absl::StreamFormat(">%s\n%s\n", pathID, path_seq);
  outStream.close();
}

void Graph::WriteDot(usize comp_id, const std::string& suffix) const {
  const DotSerializer ds(this);
  ds.WriteComponent(comp_id, suffix);
}

void Graph::WriteDot(usize comp_id, absl::Span<const PathNodeIds> flow_paths) const {
  const DotSerializer ds(this);
  ds.WriteComponent(comp_id, flow_paths);
}

void Graph::EraseNode(NodeIterator itr) {
  if (itr == nodesMap.end() || itr->second->IsMockNode()) return;

  // remove edges associated with the to be removed nodes first, then remove the node
  for (const Edge& e : *itr->second) {
    auto neighbourItr = nodesMap.find(e.GetDstID());
    if (neighbourItr == nodesMap.end()) continue;
    neighbourItr->second->EraseEdge(itr->first, ReverseEdgeKind(e.GetEdgeKind()));
  }

  nodesMap.erase(itr);
}

void Graph::EraseNode(NodeIdentifier node_id) { return EraseNode(nodesMap.find(node_id)); }

struct AlleleSpan {
  usize StartPos = 0;   // NOLINT
  usize AlleleLen = 1;  // NOLINT

  [[nodiscard]] auto GetEndPos() const noexcept -> usize { return StartPos + AlleleLen; }
};

static inline auto GetAvgBaseQual(std::string_view baseQuals, const AlleleSpan& loc) -> double {
  const auto alleleEndPos = loc.GetEndPos();
  if (baseQuals.empty() || (alleleEndPos > baseQuals.length())) return 0.0;

  auto currSum = static_cast<int>(baseQuals[loc.StartPos]);  // NOLINT
  if (loc.AlleleLen == 1) {
    return static_cast<double>(currSum);
  }

  for (usize idx = loc.StartPos + 1; idx < alleleEndPos; ++idx) {
    currSum += static_cast<int>(baseQuals[idx]);
  }

  return static_cast<double>(currSum) / static_cast<double>(baseQuals.length());
}

static inline void IncrementHpCov(HpCov& currCov, const Strand strand, const u8 haplotypeID) {
  strand == Strand::FWD ? currCov.FwdCov++ : currCov.RevCov++;
  if (haplotypeID == 1) {
    currCov.HP1++;
  } else if (haplotypeID == 2) {
    currCov.HP2++;
  } else {
    currCov.HP0++;
  }
}

void Graph::BuildVariants(absl::Span<const Transcript> transcripts, std::vector<Variant>* variants) const {
  // One Transcript has two allele hashes – one each for ref and alt allele
  //                                   SNV                    InDel
  //     HaplotypeCentered  – ----------x---------- | --------xxxxx--------
  // 'x' marks the location of the allele sequence in the haplotype
  // For each transcript, we calculate counts of exact matches for
  // any one of the three haplotype configurations for both ref & alt

  absl::flat_hash_set<usize> hapLens;                    // Set with Lengths of all haplotypes
  absl::flat_hash_map<u64, u64> hap2AlleleMap;           // Maps canonical haplotype hash to allele hash
  absl::flat_hash_map<u64, SampleHpCovs> sampleHapCovs;  // Maps allele hash to sample coverages
  absl::flat_hash_map<u64, AlleleSpan> alleleSpans;      // Maps canonical haplotype hash to allele span

  for (const Transcript& T : transcripts) {
    const auto isSNV = T.Code() == TranscriptCode::SNV;
    const auto alleles = T.GetAlleleHashes();
    const auto haplotypes = T.GetHaplotypesData();

    for (const auto& h : haplotypes) {
      const auto alSpan = isSNV ? AlleleSpan{h.hapLeftFlank, 1} : AlleleSpan{0, h.hapLen};
      alleleSpans.try_emplace(h.hapHash, alSpan);
      hapLens.insert(h.hapLen);

      const auto currAlleleHash = h.alleleKind == Allele::ALT ? alleles.AltHash : alleles.RefHash;
      sampleHapCovs.try_emplace(currAlleleHash, SampleHpCovs{});
      hap2AlleleMap.try_emplace(h.hapHash, currAlleleHash);
    }
  }

  // mateMer -> readName + sampleLabel, alleleHash
  using MateMer = std::pair<std::string, u64>;
  absl::flat_hash_set<MateMer> seenMateMers;
  const auto minBQ = static_cast<double>(params->minBaseQual);

  for (const auto& rd : sampleReads) {
    for (const usize& haplotypeLength : hapLens) {
      if (rd.Length() < haplotypeLength) continue;

      const auto readSeq = rd.SeqView();
      const auto readQual = rd.QualView();

      for (usize offset = 0; offset <= (rd.Length() - haplotypeLength); ++offset) {
        const auto merHash = Kmer::CanonicalSeqHash(absl::ClippedSubstr(readSeq, offset, haplotypeLength));
        if (!hap2AlleleMap.contains(merHash)) continue;

        const auto isTmrRd = rd.label == SampleLabel::TUMOR;
        const auto merQuals = absl::ClippedSubstr(readQual, offset, haplotypeLength);
        const auto avgAlleleQual = GetAvgBaseQual(merQuals, alleleSpans.at(merHash));

        // Skip adding to kmer count if allele length is a single base.
        // This is to reduce adding coverage from low quality bases for SNVs leading to FPs
        // Always add to kmer count for normal sample reads, so that we don't call FPs
        if (isTmrRd && avgAlleleQual < minBQ) continue;

        // Add label to key, so that keys are unique for tumor and normal samples
        const auto mmId = std::make_pair(rd.readName + ToString(rd.label), hap2AlleleMap.at(merHash));
        if (seenMateMers.find(mmId) != seenMateMers.end()) continue;

        auto& data = sampleHapCovs.at(hap2AlleleMap.at(merHash));
        IncrementHpCov(isTmrRd ? data.TmrCov : data.NmlCov, rd.strand, rd.haplotypeID);
        seenMateMers.insert(mmId);
      }
    }
  }

  usize numVariants = 0;
  for (const Transcript& T : transcripts) {
    const auto alleles = T.GetAlleleHashes();
    const auto refCovs = sampleHapCovs.at(alleles.RefHash);
    const auto altCovs = sampleHapCovs.at(alleles.AltHash);
    Variant V(T, kmerSize, VariantHpCov(refCovs.TmrCov, altCovs.TmrCov), VariantHpCov(refCovs.NmlCov, altCovs.NmlCov));

    if (V.ComputeState() == VariantState::NONE) continue;
    variants->emplace_back(std::move(V));
    numVariants++;
  }

  const auto windowId = window->ToRegionString();
  LOG_DEBUG("Built {} variants from graph for {}", numVariants, windowId);
}

auto Graph::FindRefEnd(GraphEnd k, usize comp_id, absl::Span<const NodeIdentifier> ref_mer_hashes) const
    -> Graph::RefEndResult {
  // find node in component for the first/last occurring ref kmer id.
  const auto minEndCov = static_cast<u16>(params->minAnchorCov);
  const auto numRefMers = static_cast<i64>(ref_mer_hashes.size());

  for (auto refIdx = 0; refIdx < numRefMers; ++refIdx) {
    const auto merIndex = static_cast<usize>(refIdx);
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
    if (buddyItr == nodesMap.end() || buddyItr->second == nullptr) continue;

    const auto buddysNeighbours = buddyItr->second->FindMergeableNeighbours();
    const auto areMutualBuddies = std::any_of(buddysNeighbours.cbegin(), buddysNeighbours.cend(),
                                              [&src_id](const NodeNeighbour& n) { return n.buddyId == src_id; });
    if (!areMutualBuddies) continue;

    const auto mergeDir = SourceStrand(srcNbour.edgeKind) == Strand::FWD ? BuddyPosition::FRONT : BuddyPosition::BACK;
    const auto canMergeWithBuddySeq = srcItr->second->CanMerge(*buddyItr->second, mergeDir, kmerSize);
    if (canMergeWithBuddySeq) results.emplace(srcNbour);
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
      const auto buddyNeighbourId = buddyE.GetDstID();
      if (buddyNeighbourId == src_id) continue;

      auto buddysNeighbourItr = nodesMap.find(buddyNeighbourId);
      if (buddysNeighbourItr == nodesMap.end()) continue;

      const auto srcLinkStrand = srcBuddyDiffStrands ? ReverseStrand(buddyE.GetSrcDir()) : buddyE.GetSrcDir();
      const auto resultKind = MakeEdgeKind(srcLinkStrand, buddyE.GetDstDir());

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

// NOLINTNEXTLINE
auto Graph::HasCycle(NodeIdentifier node_id, Strand direction, absl::flat_hash_set<NodeIdentifier>* touched) const
    -> bool {
  const auto itr = nodesMap.find(node_id);
  if (itr == nodesMap.end()) return false;

  touched->insert(node_id);
  for (const Edge& e : *itr->second) {
    const auto neighbourId = e.GetDstID();
    if (neighbourId == MOCK_SOURCE_ID || neighbourId == MOCK_SINK_ID || e.GetSrcDir() != direction) continue;
    if (touched->find(neighbourId) == touched->end()) return HasCycle(neighbourId, e.GetDstDir(), touched);
    touched->erase(node_id);
    return true;
  }

  touched->erase(node_id);
  return false;
}

void Graph::DisconnectEdgesTo(NodeIterator itr, const NodeContainer& nc) {
  LANCET_ASSERT(itr != nc.end());  // NOLINT

  for (const Edge& e : *itr->second) {
    auto neighbourItr = nc.find(e.GetDstID());
    if (neighbourItr == nc.end()) continue;
    neighbourItr->second->EraseEdge(itr->first);
  }
}
}  // namespace lancet2

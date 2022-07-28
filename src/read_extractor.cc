#include "lancet2/read_extractor.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/ascii.h"
#include "absl/types/span.h"
#include "lancet2/assert_macro.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#pragma clang diagnostic ignored "-Wcast-qual"
#elif defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "htslib/sam.h"

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace lancet2 {
ReadExtractor::ReadExtractor(std::shared_ptr<const CliParams> p)
    : labels({SampleLabel::TUMOR, SampleLabel::NORMAL}), params(std::move(p)) {
  auto rdrT = std::make_unique<HtsReader>(params->tumorPath, params->referencePath);
  auto rdrN = std::make_unique<HtsReader>(params->normalPath, params->referencePath);

  readers.emplace_back(std::move(rdrT));
  readers.emplace_back(std::move(rdrN));
}

auto ReadExtractor::ScanRegion(const GenomicRegion& region) -> ScanRegionResult {
  const auto tmrResult = ScanSampleRegion(0, region);
  const auto nmlResult = ScanSampleRegion(1, region);

  const auto totalNumBases = tmrResult.NumReadBases + nmlResult.NumReadBases;
  const auto avgCoverage = static_cast<double>(totalNumBases) / static_cast<double>(region.GetLength());
  const auto isActiveRegion = tmrResult.HasMutationEvidence || nmlResult.HasMutationEvidence;
  return {totalNumBases, avgCoverage, isActiveRegion};
}

auto ReadExtractor::ExtractReads(const GenomicRegion& region, double sampleFraction) -> ReadInfoList {
  std::vector<ReadInfo> finalReads;
  sampler.SetFraction(sampleFraction);

  const auto tmrMateInfo = FetchReads(0, region, &finalReads);
  if (params->extractReadPairs && !tmrMateInfo.empty()) {
    FetchPairs(0, tmrMateInfo, &finalReads);
  }

  const auto nmlMateInfo = FetchReads(1, region, &finalReads);
  if (params->extractReadPairs && !nmlMateInfo.empty()) {
    FetchPairs(1, nmlMateInfo, &finalReads);
  }

  return finalReads;
}

auto ReadExtractor::FetchReads(usize sampleIdx, const GenomicRegion& region, ReadInfoList* result) -> MateInfoMap {
  HtsReader* rdr = readers.at(sampleIdx).get();
  const auto label = labels.at(sampleIdx);

  const auto jumpStatus = rdr->JumpToRegion(region);
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  // readName -> mateRegion
  absl::flat_hash_map<std::string, GenomicRegion> mateName2Region;
  HtsAlignment aln;

  while (rdr->GetNextAlignment(&aln, {"XT", "XA", "AS", "XS"}) == HtsReader::IteratorState::VALID) {
    if (!sampler.ShouldSample() || !PassesFilters(aln, *params, label)) continue;

    auto rdInfo = aln.BuildReadInfo(label, params->trimBelowQual, params->maxKmerSize);
    if (rdInfo.IsEmpty()) continue;

    if (params->extractReadPairs && !aln.IsMateUnmapped()) {
      const auto itr = mateName2Region.find(aln.GetReadName());
      if (itr == mateName2Region.end()) {
        mateName2Region.emplace(aln.GetReadName(), aln.GetMateRegion());
      } else {
        mateName2Region.erase(itr);
      }
    }

    result->emplace_back(std::move(rdInfo));
  }

  return mateName2Region;
}

void ReadExtractor::FetchPairs(usize sampleIdx, const MateInfoMap& mate_info, ReadInfoList* result) {
  HtsReader* rdr = readers.at(sampleIdx).get();
  const auto label = labels.at(sampleIdx);
  // sort mate positions map in co-ordinate sorted order
  std::vector<GenomicRegion> mateRegions;
  mateRegions.reserve(mate_info.size());
  std::for_each(mate_info.cbegin(), mate_info.cend(),
                [&mateRegions](const auto& kv) { mateRegions.emplace_back(kv.second); });
  std::sort(mateRegions.begin(), mateRegions.end(), [&rdr](const GenomicRegion& gr1, const GenomicRegion& gr2) -> bool {
    if (gr1.GetChromName() != gr2.GetChromName()) {
      return rdr->GetContigIndex(gr1.GetChromName()) < rdr->GetContigIndex(gr2.GetChromName());
    }

    if (gr1.GetStartPos1() != gr2.GetStartPos1()) return gr1.GetStartPos1() < gr2.GetStartPos1();
    return gr1.GetEndPos1() < gr2.GetEndPos1();
  });

  const auto jumpStatus = rdr->SetBatchRegions(absl::MakeConstSpan(mateRegions));
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  HtsAlignment aln;
  while (rdr->GetNextAlignment(&aln, {}) == HtsReader::IteratorState::VALID) {
    if (!PassesFilters(aln, *params, label) || !mate_info.contains(aln.GetReadName())) continue;

    auto rdInfo = aln.BuildReadInfo(label, params->trimBelowQual, params->maxKmerSize);
    if (rdInfo.IsEmpty()) continue;

    result->emplace_back(std::move(rdInfo));
  }
}

auto ReadExtractor::PassesFilters(const HtsAlignment& aln, const CliParams& params, SampleLabel label) -> bool {
  if (label == SampleLabel::NORMAL) return true;

  if (aln.IsQcFailed() || aln.IsDuplicate()) return false;
  if (params.skipSecondary && aln.IsSecondary()) return false;

  // XT type: Unique/Repeat/N/Mate-sw
  // XT:A:M (one-mate recovered) means that one of the pairs is uniquely mapped and the other isn't
  // Heng Li: If the read itself is a repeat and can't be mapped without relying on its mate, you
  // see "XT:Z:R". Nonetheless, the mapping quality is not necessarily zero. When its mate can be
  // mapped unambiguously, the read can still be mapped confidently and thus assigned a high mapQ.
  // MapQ is computed for the read pair. XT is determined from a single read.
  if (aln.HasTag("XT")) return false;

  // XA -- BWA (Illumina): alternative hits; format: (chr,pos,CIGAR,NM;)
  if (aln.HasTag("XA")) return false;

  if (aln.GetMappingQual() < params.minReadMappingQual) return false;

  // AS: Alignment score
  // XS: Suboptimal alignment score
  if (aln.HasTag("AS") && aln.HasTag("XS")) {
    const auto AS = bam_aux2i(aln.GetTagData("AS").value());
    const auto XS = bam_aux2i(aln.GetTagData("XS").value());
    return std::abs(AS - XS) >= params.minReadAsXsDiff;
  }

  return true;
}

auto ReadExtractor::ScanSampleRegion(usize sampleIdx, const GenomicRegion& region) -> ScanRegionResult {
  HtsReader* rdr = readers.at(sampleIdx).get();
  const auto label = labels.at(sampleIdx);
  const auto jumpStatus = rdr->JumpToRegion(region);
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  using u32 = u32;
  std::vector<u32> genomePositions;  // softclip genome positions for single alignment
  std::map<u32, u32> mismatches;     // genome position -> number of mismatches at position
  std::map<u32, u32> insertions;     // genome position -> number of insertions at position
  std::map<u32, u32> deletions;      // genome position -> number of deletions at position
  std::map<u32, u32> softclips;      // genome position -> number of softclips at position

  // lambda function to increment the counter for `genomePositions`
  const auto incrementGPos = [](std::map<u32, u32>& m, u32 genome_pos) {
    auto&& itr = m.find(genome_pos);
    if (itr != m.end()) {
      ++itr->second;
      return;
    }
    m.emplace(std::make_pair(genome_pos, 1));
  };

  // lambda function to check if map `m` has evidence of mutation >= minAltCntTumor
  const auto hasEvidence = [](const std::map<u32, u32>& m, u32 min_alt_cnt_tmr) {
    return std::any_of(m.cbegin(), m.cend(), [&min_alt_cnt_tmr](const auto& p) { return p.second >= min_alt_cnt_tmr; });
  };

  HtsAlignment aln;
  bool isActiveRegion = false;
  u64 numReadBases = 0;

  while (rdr->GetNextAlignment(&aln, {"MD"}) == HtsReader::IteratorState::VALID) {
    if (!aln.IsUnmapped() && !aln.IsDuplicate()) numReadBases += aln.GetLength();

    // skip processing further if already an active region or if it doesn't pass filters
    if (isActiveRegion || !PassesFilters(aln, *params, label)) continue;

    if (aln.HasTag("MD")) {
      FillMD(bam_aux2Z(aln.GetTagData("MD").value()), aln.GetReadQuality(), aln.GetStartPos0(), params->minBaseQual,
             &mismatches);
    }

    const auto cigarUnits = aln.GetCigarData();
    auto currGPos = static_cast<u32>(aln.GetStartPos0());
    bool hasSoftClip = false;

    for (const auto& cigUnit : cigarUnits) {
      if (cigUnit.ConsumesReference()) currGPos += cigUnit.Length;
      switch (cigUnit.Operation) {
        case CigarOp::INSERTION:
          incrementGPos(insertions, currGPos);
          LANCET_ASSERT(!insertions.empty());  // NOLINT
          break;
        case CigarOp::DELETION:
          incrementGPos(deletions, currGPos);
          LANCET_ASSERT(!deletions.empty());  // NOLINT
          break;
        case CigarOp::SEQUENCE_MISMATCH:
          incrementGPos(mismatches, currGPos);
          LANCET_ASSERT(!mismatches.empty());  // NOLINT
          break;
        case CigarOp::SOFT_CLIP:
          hasSoftClip = true;
          break;
        default:
          break;
      }
    }

    genomePositions.clear();
    if (hasSoftClip && aln.GetSoftClips(nullptr, nullptr, &genomePositions, false)) {
      for (const auto gpos : genomePositions) incrementGPos(softclips, gpos);
      LANCET_ASSERT(!softclips.empty());  // NOLINT
    }

    // if evidence of SNV/insertion/deletion/softclip is found in window, mark region as active
    isActiveRegion = hasEvidence(mismatches, params->minTmrAltCnt) || hasEvidence(insertions, params->minTmrAltCnt) ||
                     hasEvidence(deletions, params->minTmrAltCnt) || hasEvidence(softclips, params->minTmrAltCnt);
  }

  auto avgCoverage = static_cast<double>(numReadBases) / static_cast<double>(region.GetLength());
  isActiveRegion = hasEvidence(mismatches, params->minTmrAltCnt) || hasEvidence(insertions, params->minTmrAltCnt) ||
                   hasEvidence(deletions, params->minTmrAltCnt) || hasEvidence(softclips, params->minTmrAltCnt);

  return ScanRegionResult{numReadBases, avgCoverage, isActiveRegion};
}

void ReadExtractor::FillMD(std::string_view md, std::string_view quals, i64 aln_start, u32 min_bq,
                           std::map<u32, u32>* result) {
  if (aln_start < 0) return;
  auto genomePos = static_cast<u32>(aln_start);
  std::string token;

  for (const auto& c : md) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(c))) {
      token += c;
      continue;
    }

    const auto step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genomePos += static_cast<u32>(step);
    token.clear();

    const auto basePos = static_cast<usize>(genomePos - aln_start);
    if (static_cast<int>(quals[basePos]) < static_cast<int>(min_bq)) continue;

    const auto base = absl::ascii_toupper(static_cast<unsigned char>(c));
    if (base == 'A' || base == 'C' || base == 'T' || base == 'G') {
      const auto& itr = result->find(genomePos);
      if (itr != result->end()) {
        itr->second++;
      } else {
        result->emplace(genomePos, 1);
      }
    }
  }
}
}  // namespace lancet2

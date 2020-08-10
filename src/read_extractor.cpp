#include "lancet/read_extractor.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/ascii.h"
#include "absl/types/span.h"
#include "htslib/sam.h"
#include "lancet/assert_macro.h"
#include "lancet/fractional_sampler.h"

namespace lancet {
ReadExtractor::ReadExtractor(std::shared_ptr<const CliParams> p)
    : tmrRdr(p->tumorPath, p->referencePath), nmlRdr(p->normalPath, p->referencePath) {
  params = std::move(p);
}

void ReadExtractor::SetTargetRegion(const GenomicRegion& region) {
  targetRegion = region;
  const auto tmrResult = EvaluateRegion(&tmrRdr, targetRegion, *params);
  const auto nmlResult = EvaluateRegion(&nmlRdr, targetRegion, *params);

  isActiveRegion = tmrResult.isActiveRegion || nmlResult.isActiveRegion;
  avgCoverage = (tmrResult.coverage + nmlResult.coverage) / 2.0F;
}

auto ReadExtractor::Extract() -> ReadInfoList {
  std::vector<ReadInfo> finalReads;
  const auto tmrMateInfo = ExtractReads(&tmrRdr, &finalReads, SampleLabel::TUMOR);
  if (params->extractReadPairs && !tmrMateInfo.empty()) {
    ExtractPairs(&tmrRdr, tmrMateInfo, &finalReads, SampleLabel::TUMOR);
  }

  const auto nmlMateInfo = ExtractReads(&nmlRdr, &finalReads, SampleLabel::NORMAL);
  if (params->extractReadPairs && !nmlMateInfo.empty()) {
    ExtractPairs(&nmlRdr, nmlMateInfo, &finalReads, SampleLabel::NORMAL);
  }

  return finalReads;
}

auto ReadExtractor::ExtractReads(HtsReader* rdr, ReadInfoList* result, SampleLabel label)
    -> absl::flat_hash_map<std::string, GenomicRegion> {
  const auto jumpStatus = rdr->SetRegion(targetRegion);
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  const auto fractionToSample = avgCoverage > params->maxWindowCov ? (avgCoverage / params->maxWindowCov) : 1.0;
  FractionalSampler sampler(fractionToSample);

  // readName -> mateRegion
  absl::flat_hash_map<std::string, GenomicRegion> mateName2Region;
  HtsAlignment aln;

  while (rdr->NextAlignment(&aln, {"XT", "XA", "AS", "XS"}) == HtsReader::IteratorState::VALID) {
    if (!sampler.ShouldSample() || !PassesFilters(aln, *params)) continue;
    if (!params->useOverlapReads && !aln.IsWithinRegion(targetRegion)) continue;
    if (label == SampleLabel::TUMOR && !PassesTmrFilters(aln, *params)) continue;

    auto rdInfo = aln.BuildReadInfo(label, params->trimBelowQual, params->maxKmerSize);
    if (rdInfo.IsEmpty()) continue;

    if (params->extractReadPairs && !aln.IsMateUnmapped()) {
      const auto itr = mateName2Region.find(aln.ReadName());
      if (itr == mateName2Region.end()) {
        mateName2Region.emplace(aln.ReadName(), aln.MateRegion());
      } else {
        mateName2Region.erase(itr);
      }
    }

    result->emplace_back(std::move(rdInfo));
  }

  return mateName2Region;
}

void ReadExtractor::ExtractPairs(HtsReader* rdr, const absl::flat_hash_map<std::string, GenomicRegion>& mate_info,
                                 ReadInfoList* result, SampleLabel label) {
  // sort mate positions map in co-ordinate sorted order
  std::vector<GenomicRegion> mateRegions;
  mateRegions.reserve(mate_info.size());
  std::for_each(mate_info.cbegin(), mate_info.cend(),
                [&mateRegions](const auto& kv) { mateRegions.emplace_back(kv.second); });
  std::sort(mateRegions.begin(), mateRegions.end(), [&rdr](const GenomicRegion& gr1, const GenomicRegion& gr2) -> bool {
    if (gr1.Chromosome() != gr2.Chromosome()) return rdr->ContigID(gr1.Chromosome()) < rdr->ContigID(gr2.Chromosome());
    if (gr1.StartPosition1() != gr2.StartPosition1()) return gr1.StartPosition1() < gr2.StartPosition1();
    return gr1.EndPosition1() < gr2.EndPosition1();
  });

  const auto jumpStatus = rdr->SetRegions(absl::MakeConstSpan(mateRegions));
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  HtsAlignment aln;
  while (rdr->NextAlignment(&aln, {}) == HtsReader::IteratorState::VALID) {
    if (!PassesFilters(aln, *params) || !mate_info.contains(aln.ReadName())) continue;
    result->emplace_back(aln.BuildReadInfo(label, params->trimBelowQual, params->maxKmerSize));
  }
}

auto ReadExtractor::PassesFilters(const HtsAlignment& aln, const CliParams& params) -> bool {
  if (aln.IsQcFailed() || aln.IsDuplicate()) return false;
  return !(params.skipSecondary && aln.IsSecondary());
}

auto ReadExtractor::PassesTmrFilters(const HtsAlignment& aln, const CliParams& params) -> bool {
  // XT type: Unique/Repeat/N/Mate-sw
  // XT:A:M (one-mate recovered) means that one of the pairs is uniquely mapped and the other isn't
  // Heng Li: If the read itself is a repeat and can't be mapped without relying on its mate, you
  // see "XT:Z:R". Nonetheless, the mapping quality is not necessarily zero. When its mate can be
  // mapped unambiguously, the read can still be mapped confidently and thus assigned a high mapQ.
  // MapQ is computed for the read pair. XT is determined from a single read.
  if (aln.HasTag("XT")) return false;

  // XA -- BWA (Illumina): alternative hits; format: (chr,pos,CIGAR,NM;)
  if (aln.HasTag("XA")) return false;

  if (aln.MappingQuality() < params.minReadMappingQual) return false;

  // AS: Alignment score
  // XS: Suboptimal alignment score
  if (aln.HasTag("AS") && aln.HasTag("XS")) {
    const auto AS = bam_aux2i(aln.TagData("AS").ValueOrDie());
    const auto XS = bam_aux2i(aln.TagData("XS").ValueOrDie());
    return std::abs(AS - XS) >= params.minReadAsXsDiff;
  }

  return true;
}

auto ReadExtractor::EvaluateRegion(HtsReader* rdr, const GenomicRegion& region, const CliParams& params) -> EvalResult {
  const auto jumpStatus = rdr->SetRegion(region);
  if (!jumpStatus.ok()) throw std::runtime_error(jumpStatus.ToString());

  using u32 = std::uint32_t;
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
  std::uint64_t numReadBases = 0;

  while (rdr->NextAlignment(&aln, {"MD"}) == HtsReader::IteratorState::VALID) {
    if (!aln.IsUnmapped() && !aln.IsDuplicate()) numReadBases += aln.Length();

    // skip processing further if already an active region or if doesn't pass filters
    if (isActiveRegion || !PassesFilters(aln, params)) continue;
    if (!params.useOverlapReads && !aln.IsWithinRegion(region)) continue;

    if (aln.HasTag("MD")) {
      FillMDMismatches(bam_aux2Z(aln.TagData("MD").ValueOrDie()), aln.ReadQuality(), aln.StartPosition0(),
                       params.minBaseQual, &mismatches);
    }

    const auto cigarUnits = aln.CigarData();
    auto currGPos = static_cast<u32>(aln.StartPosition0());
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
    if (hasSoftClip && aln.SoftClips(nullptr, nullptr, &genomePositions, false)) {
      for (const auto gpos : genomePositions) incrementGPos(softclips, gpos);
      LANCET_ASSERT(!softclips.empty());  // NOLINT
    }

    // if evidence of SNV/insertion/deletion/softclip is found in window, mark region as active
    isActiveRegion = hasEvidence(mismatches, params.minTmrAltCnt) || hasEvidence(insertions, params.minTmrAltCnt) ||
                     hasEvidence(deletions, params.minTmrAltCnt) || hasEvidence(softclips, params.minTmrAltCnt);
  }

  auto avgCoverage = static_cast<double>(numReadBases) / static_cast<double>(region.Length());
  isActiveRegion = hasEvidence(mismatches, params.minTmrAltCnt) || hasEvidence(insertions, params.minTmrAltCnt) ||
                   hasEvidence(deletions, params.minTmrAltCnt) || hasEvidence(softclips, params.minTmrAltCnt);

  return EvalResult{avgCoverage, isActiveRegion};
}

void ReadExtractor::FillMDMismatches(std::string_view md, std::string_view quals, std::int64_t aln_start,
                                     std::uint32_t min_bq, std::map<std::uint32_t, std::uint32_t>* result) {
  if (aln_start < 0) return;
  auto genomePos = static_cast<std::uint32_t>(aln_start);
  std::string token;

  for (const auto& c : md) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(c))) {
      token += c;
      continue;
    }

    const auto step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genomePos += static_cast<std::uint32_t>(step);
    token.clear();

    const auto basePos = static_cast<std::size_t>(genomePos - aln_start);
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
}  // namespace lancet

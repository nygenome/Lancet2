#include "lancet/core/read_collector.h"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string>
#include <utility>

#include "absl/strings/ascii.h"
#include "lancet/base/assert.h"
#include "lancet/base/online_stats.h"
#include "spdlog/fmt/fmt.h"

namespace {

inline void ParseMd(std::string_view md_val, absl::Span<const u8> quals, const i64 start, std::map<u32, u32>* result) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (start < 0) return;

  auto genome_pos = static_cast<u32>(start);
  std::string token;

  for (const auto& character : md_val) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(character))) {
      token += character;
      continue;
    }

    const auto step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genome_pos += static_cast<u32>(step);
    token.clear();

    const auto base_pos = static_cast<usize>(genome_pos - start);
    static constexpr u8 min_bq = 20;
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (quals.at(base_pos) < min_bq) continue;

    const auto base = absl::ascii_toupper(static_cast<unsigned char>(character));
    if (base == 'A' || base == 'C' || base == 'T' || base == 'G') {
      const auto& itr = result->find(genome_pos);
      if (itr != result->end()) {
        itr->second++;
      } else {
        result->emplace(genome_pos, 1);
      }
    }
  }
}

[[nodiscard]] inline auto HasOnlyMatches(absl::Span<const lancet::hts::CigarUnit> cigar) -> bool {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (cigar.size() != 1) return false;

  using lancet::hts::CigarOp;
  const auto cigop = cigar[0].Operation();
  return (cigop == CigarOp::ALIGNMENT_MATCH || cigop == CigarOp::SEQUENCE_MATCH);
}

}  // namespace

namespace lancet::core {

ReadCollector::ReadCollector(Params params) : mParams(std::move(params)) {
  using hts::Extractor;
  using hts::Alignment::Fields::AUX_RGAUX;
  using hts::Alignment::Fields::CIGAR_SEQ_QUAL;
  static const auto aln_tags = std::vector<std::string>{"AS", "XS", "XT", "XA", "SA"};
  const auto needs_pairs = mParams.mExtractReadPairs;

  mSampleList = MakeSampleList(mParams);
  for (const auto& sinfo : mSampleList) {
    const auto needs_tags = (sinfo.TagKind() == cbdg::Label::TUMOR && !(mParams.mNoFilterRds)) || needs_pairs;
    const auto fields = needs_tags ? AUX_RGAUX : CIGAR_SEQ_QUAL;
    const auto tags = needs_tags ? absl::MakeConstSpan(aln_tags) : absl::Span<const std::string>();
    auto extractor = std::make_unique<Extractor>(sinfo.Path(), mParams.mRefPath, fields, tags, mParams.mNoCtgCheck);
    mExtractors.emplace(sinfo, std::move(extractor));
  }
}

auto ReadCollector::CollectRegionResult(const Region& region) -> Result {
  const auto max_sample_cov = mParams.mMaxWinCov / static_cast<f64>(mSampleList.size());

  // custom hash and equal, so we don't add reads with same seq multiple times
  std::vector<Read> sample_reads;
  absl::flat_hash_map<std::string, hts::Alignment::MateInfo> expected_mate_regions;

  for (auto& sinfo : mSampleList) {
    u64 num_reads = 0;
    u64 num_bases = 0;
    expected_mate_regions.clear();

    const AlnAndRefPaths aln_refs{sinfo.Path(), mParams.mRefPath};
    const auto sample_name = std::string(sinfo.SampleName());
    const auto is_tumor_sample = sinfo.TagKind() == cbdg::Label::TUMOR;
    const auto sample_cov = EstimateCoverage(sinfo, region);

    const auto pct_to_sample = sample_cov > max_sample_cov ? ((max_sample_cov * 100.0) / sample_cov) : 100.0;
    mDownsampler.SetPercentToSample(pct_to_sample);

    auto& extractor = mExtractors.at(sinfo);
    extractor->SetRegionToExtract(region.ToSamtoolsRegion());

    for (const auto& aln : *extractor) {
      const auto bflag = aln.Flag();
      // NOLINTBEGIN(readability-braces-around-statements)
      if (bflag.IsDuplicate() || bflag.IsQcFail() || bflag.IsSecondary()) continue;
      if (is_tumor_sample && !mParams.mNoFilterRds && FailsFilter(aln)) continue;
      if (!mDownsampler.ShouldSample()) continue;
      // NOLINTEND(readability-braces-around-statements)

      num_reads += 1;
      num_bases += aln.Length();
      sample_reads.emplace_back(aln, sample_name, sinfo.TagKind());

      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!mParams.mExtractReadPairs) continue;

      // First check if we already saw both mates in the same window
      const auto mate_itr = expected_mate_regions.find(aln.QnameView());
      if (mate_itr != expected_mate_regions.end()) {
        expected_mate_regions.erase(mate_itr);
        continue;
      }

      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (bflag.IsMateUnmapped()) continue;

      const auto has_split_aln = aln.HasTag("SA");
      const auto curr_insert = std::abs(aln.InsertSize());
      const auto abnormal_insert = curr_insert < sinfo.mMinExpectedInsert || curr_insert > sinfo.mMaxExpectedInsert;
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!abnormal_insert && !has_split_aln && bflag.IsMappedProperPair()) continue;

      auto [itr, newly_emplaced] = expected_mate_regions.try_emplace(aln.QnameView(), aln.MateLocation());
      // If not newly emplaced, then we already read both pairs
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!newly_emplaced) expected_mate_regions.erase(itr);
    }

    if (!mParams.mExtractReadPairs || expected_mate_regions.empty()) {
      sinfo.SetNumReads(num_reads);
      sinfo.SetNumBases(num_bases);
      sinfo.CalculateMeanCov(region.Length());
      continue;
    }

    std::vector<std::string> mate_region_specs;
    mate_region_specs.reserve(expected_mate_regions.size());
    std::ranges::transform(BuildSortedMateInfos(expected_mate_regions), std::back_inserter(mate_region_specs),
                           [&extractor](const hts::Alignment::MateInfo& item) -> std::string {
                             const auto mate_chrom = extractor->ChromName(item.mChromIndex);
                             const auto mate_pos1 = item.mMateStartPos0 + 1;
                             const auto colon_in_mate_chrom = mate_chrom.find(':') != std::string::npos;
                             return colon_in_mate_chrom ? fmt::format("{{{}}}:{}-{}", mate_chrom, mate_pos1, mate_pos1)
                                                        : fmt::format("{}:{}-{}", mate_chrom, mate_pos1, mate_pos1);
                           });

    for (const auto& mr_spec : mate_region_specs) {
      extractor->SetRegionToExtract(mr_spec);

      for (const auto& aln : *extractor) {
        const auto itr = expected_mate_regions.find(aln.QnameView());
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (itr == expected_mate_regions.end()) continue;

        num_reads += 1;
        num_bases += aln.Length();
        sample_reads.emplace_back(aln, sample_name, sinfo.TagKind());

        // We break out of current mate region when we find the expected mate
        expected_mate_regions.erase(itr);
        break;
      }
    }

    sinfo.SetNumReads(num_reads);
    sinfo.SetNumBases(num_bases);
    sinfo.CalculateMeanCov(region.Length());
  }

  std::ranges::sort(sample_reads, [](const Read& lhs, const Read& rhs) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (lhs.TagKind() != rhs.TagKind()) return static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind());
    if (lhs.SampleName() != rhs.SampleName()) return lhs.SampleName() < rhs.SampleName();
    if (lhs.QnameView() != rhs.QnameView()) return lhs.QnameView() < rhs.QnameView();
    if (lhs.ChromIndex() != rhs.ChromIndex()) return lhs.ChromIndex() < rhs.ChromIndex();
    return lhs.StartPos0() < rhs.StartPos0();
    // NOLINTEND(readability-braces-around-statements)
  });

  return {.mSampleReads = std::move(sample_reads), .mSampleList = mSampleList};
}

auto ReadCollector::EstimateInsertRange(const AlnAndRefPaths& paths) -> InsertRange {
  static constexpr usize MAX_READS_TO_SAMPLE = 100000;
  static constexpr usize MIN_REQUIRED_MAPPING_QUALITY = 30;

  OnlineStats stats;
  InsertRange normal_range{0, 0};
  const auto [bam_cram_path, ref_path] = paths;
  hts::Extractor extractor(bam_cram_path, ref_path, hts::Alignment::Fields::CIGAR_SEQ_QUAL, {}, true);

  for (const auto& aln : extractor) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (stats.Count() >= MAX_READS_TO_SAMPLE) break;

    const auto bitflag = aln.Flag();
    if (bitflag.IsSecondary() || bitflag.IsSupplementary() || bitflag.IsQcFail() || bitflag.IsDuplicate() ||
        aln.ChromIndex() != aln.MateChromIndex() || aln.MapQual() < MIN_REQUIRED_MAPPING_QUALITY) {
    }

    const auto cigar = aln.CigarData();
    if (bitflag.IsPairedInSequencing() && bitflag.IsMappedProperPair() && HasOnlyMatches(cigar)) {
      stats.Add(std::abs(aln.InsertSize()));
    }
  }

  const auto mean_iss = stats.Mean();
  const auto stddev_iss = stats.StdDev();

  static constexpr f64 SIGMA = 2.0;
  normal_range[0] = static_cast<i64>(mean_iss - (SIGMA * stddev_iss));
  normal_range[1] = static_cast<i64>(mean_iss + (SIGMA * stddev_iss));

  return normal_range;
}

auto ReadCollector::IsActiveRegion(const Params& params, const Region& region) -> bool {
  std::vector<u32> genome_positions;  // softclip genome positions for single alignment
  std::map<u32, u32> mismatches;      // genome position -> number of mismatches at position
  std::map<u32, u32> insertions;      // genome position -> number of insertions at position
  std::map<u32, u32> deletions;       // genome position -> number of deletions at position
  std::map<u32, u32> softclips;       // genome position -> number of softclips at position

  for (const auto& sinfo : MakeSampleList(params)) {
    genome_positions.clear();
    mismatches.clear();
    insertions.clear();
    deletions.clear();
    softclips.clear();

    constexpr auto fields = hts::Alignment::Fields::AUX_RGAUX;
    hts::Extractor extractor(sinfo.Path(), params.mRefPath, fields, {"MD"}, params.mNoCtgCheck);
    extractor.SetRegionToExtract(region.ToSamtoolsRegion());

    for (const auto& aln : extractor) {
      const auto bflag = aln.Flag();
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (bflag.IsDuplicate() || bflag.IsQcFail() || bflag.IsSecondary() || bflag.IsUnmapped()) continue;

      if (aln.HasTag("MD")) {
        ParseMd(aln.GetTag<std::string_view>("MD").value(), aln.QualView(), aln.StartPos0(), &mismatches);
      }

      const auto cigar_units = aln.CigarData();
      auto curr_genome_pos = static_cast<u32>(aln.StartPos0());
      bool has_soft_clip = false;

      // lambda function to increment the counter for `genome_positions`
      static const auto increment_genome_pos = [](std::map<u32, u32>& map, u32 genome_pos) {
        auto&& itr = map.find(genome_pos);
        if (itr != map.end()) {
          ++itr->second;
          return;
        }
        map.emplace(genome_pos, 1);
      };

      for (const auto& cig_unit : cigar_units) {
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (cig_unit.ConsumesReference()) curr_genome_pos += cig_unit.Length();
        switch (cig_unit.Operation()) {
          case hts::CigarOp::INSERTION:
            increment_genome_pos(insertions, curr_genome_pos);
            LANCET_ASSERT(!insertions.empty())
            break;
          case hts::CigarOp::DELETION:
            increment_genome_pos(deletions, curr_genome_pos);
            LANCET_ASSERT(!deletions.empty())
            break;
          case hts::CigarOp::SEQUENCE_MISMATCH:
            increment_genome_pos(mismatches, curr_genome_pos);
            LANCET_ASSERT(!mismatches.empty())
            break;
          case hts::CigarOp::SOFT_CLIP:
            has_soft_clip = true;
            break;
          default:
            break;
        }
      }

      genome_positions.clear();
      if (has_soft_clip && aln.GetSoftClips(nullptr, nullptr, &genome_positions, false)) {
        // NOLINTNEXTLINE(readability-braces-around-statements)
        for (const auto gpos : genome_positions) increment_genome_pos(softclips, gpos);
        LANCET_ASSERT(!softclips.empty())
      }

      // lambda function to check if `map` has evidence of mutation >= 2 reads
      static const auto found_mutation_evidence = [](const std::map<u32, u32>& map) {
        return std::any_of(map.cbegin(), map.cend(), [](const auto& pair) { return pair.second >= 2; });
      };

      // if evidence of SNV/insertion/deletion/softclip is found in window, mark region as active
      if (found_mutation_evidence(mismatches) || found_mutation_evidence(insertions) ||
          found_mutation_evidence(deletions) || found_mutation_evidence(softclips)) {
        return true;
      }
    }
  }

  return false;
}

auto ReadCollector::BuildSampleNameList(const Params& params) -> std::vector<std::string> {
  const auto sinfo_list = MakeSampleList(params);
  std::vector<std::string> results;
  results.reserve(sinfo_list.size());
  std::ranges::transform(sinfo_list, std::back_inserter(results),
                         [](const SampleInfo& item) -> std::string { return std::string(item.SampleName()); });
  return results;
}

auto ReadCollector::EstimateCoverage(const SampleInfo& sinfo, const Region& region) const -> f64 {
  static const auto filter_tags = std::vector<std::string>{"AS", "XS", "XT", "XA"};

  const auto need_pairs = mParams.mExtractReadPairs;
  const auto need_filt = sinfo.TagKind() == cbdg::Label::NORMAL || !(mParams.mNoFilterRds);
  const auto fields = need_filt ? hts::Alignment::Fields::AUX_RGAUX : hts::Alignment::Fields::CORE_QNAME;
  const auto fill_tags = need_filt ? filter_tags : std::vector<std::string>{};

  hts::Extractor extractor(sinfo.Path(), mParams.mRefPath, fields, fill_tags, true);
  extractor.SetRegionToExtract(region.ToSamtoolsRegion());

  static const auto base_summer = [&need_filt, &need_pairs, &region](const u64 sum, const hts::Alignment& aln) -> u64 {
    const auto bflag = aln.Flag();
    // NOLINTBEGIN(readability-braces-around-statements)
    if (bflag.IsDuplicate() || bflag.IsQcFail() || bflag.IsSecondary() || (need_filt && FailsFilter(aln))) return sum;
    if (!need_pairs || aln.MateOverlapsRegion(region)) return sum + aln.Length();
    // NOLINTEND(readability-braces-around-statements)
    return sum + aln.Length() + aln.Length();
  };

  const u64 num_bases = std::accumulate(extractor.begin(), extractor.end(), 0, base_summer);
  return static_cast<f64>(num_bases) / static_cast<f64>(region.Length());
}

auto ReadCollector::FailsFilter(const hts::Alignment& aln) -> bool {
  static constexpr u8 DEFAULT_MIN_READ_MAPPING_QUALITY = 20;
  static constexpr i64 DEFAULT_MIN_READ_AS_XS_DIFF = 5;

  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (aln.MapQual() < DEFAULT_MIN_READ_MAPPING_QUALITY) return true;

  // AS: Alignment score
  // XS: Suboptimal alignment score
  if (aln.HasTag("AS") && aln.HasTag("XS")) {
    const auto as_tag = aln.GetTag<i64>("AS").value();
    const auto xs_tag = aln.GetTag<i64>("XS").value();
    if (std::abs(as_tag - xs_tag) < DEFAULT_MIN_READ_AS_XS_DIFF) {
      return true;
    }
  }

  // XT type: Unique/Repeat/N/Mate-sw
  // XT:A:M (one-mate recovered) means that one of the pairs is uniquely mapped and the other isn't
  // Heng Li: If the read itself is a repeat and can't be mapped without relying on its mate, you
  // see "XT:Z:R". Nonetheless, the mapping quality is not necessarily zero. When its mate can be
  // mapped unambiguously, the read can still be mapped confidently and thus assigned a high mapQ.
  // MapQ is computed for the read pair. XT is determined from a single read.
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (aln.HasTag("XT")) return true;

  // XA -- BWA (Illumina): alternative hits; format: (chr,pos,CIGAR,NM;)
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (aln.HasTag("XA")) return true;

  return false;
}

auto ReadCollector::MakeSampleList(const Params& params) -> std::vector<SampleInfo> {
  std::vector<SampleInfo> results;
  results.reserve(params.mNormals.size() + params.mTumors.size());
  const hts::Reference ref(params.mRefPath);

  // Add normal samples
  std::ranges::transform(params.mNormals, std::back_inserter(results), [&ref](const BamCramWithInsert& item) {
    const hts::Extractor extractor(item.first, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
    auto result = SampleInfo(extractor.SampleName(), item.first, cbdg::Label::NORMAL);
    result.mMinExpectedInsert = item.second[0];
    result.mMaxExpectedInsert = item.second[1];
    return result;
  });

  // Add tumor samples
  std::ranges::transform(params.mTumors, std::back_inserter(results), [&ref](const BamCramWithInsert& item) {
    const hts::Extractor extractor(item.first, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
    auto result = SampleInfo(extractor.SampleName(), item.first, cbdg::Label::TUMOR);
    result.mMinExpectedInsert = item.second[0];
    result.mMaxExpectedInsert = item.second[1];
    return result;
  });

  std::ranges::sort(results, std::less<SampleInfo>{});
  return results;
}

auto ReadCollector::BuildSortedMateInfos(const MateRegionsMap& data) -> std::vector<hts::Alignment::MateInfo> {
  std::vector<hts::Alignment::MateInfo> results;
  results.reserve(data.size());

  using MateNameAndRegion = std::pair<std::string, hts::Alignment::MateInfo>;
  std::ranges::transform(data, std::back_inserter(results), [](const MateNameAndRegion& item) { return item.second; });

  std::ranges::sort(results, [](const hts::Alignment::MateInfo& lhs, const hts::Alignment::MateInfo& rhs) -> bool {
    return (lhs.mChromIndex != rhs.mChromIndex) ? lhs.mChromIndex < rhs.mChromIndex
                                                : lhs.mMateStartPos0 < rhs.mMateStartPos0;
  });

  return results;
}

}  // namespace lancet::core

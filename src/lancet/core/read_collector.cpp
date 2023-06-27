#include "lancet/core/read_collector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <string>
#include <utility>

#include "absl/container/btree_map.h"
#include "absl/strings/ascii.h"
#include "lancet/base/assert.h"
#include "lancet/base/compute_stats.h"
#include "spdlog/fmt/fmt.h"

using CountMap = absl::btree_map<u32, u32>;

namespace {

inline void ParseMd(std::string_view md_val, absl::Span<const u8> quals, const i64 start, CountMap* result) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (start < 0) return;

  std::string token;
  token.reserve(md_val.length());
  auto genome_pos = static_cast<u32>(start);

  for (const auto& character : md_val) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(character))) {
      token += character;
      continue;
    }

    const auto step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genome_pos += static_cast<u32>(step);
    token.clear();

    const auto base_pos = static_cast<usize>(genome_pos - start);
    static constexpr u8 MIN_BASE_QUAL = 20;
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (quals.at(base_pos) < MIN_BASE_QUAL) continue;

    const auto base = absl::ascii_toupper(static_cast<unsigned char>(character));
    if (base == 'A' || base == 'C' || base == 'T' || base == 'G') {
      auto [itr, newly_added] = result->try_emplace(genome_pos, 1);
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!newly_added) itr->second++;
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

// NOLINTNEXTLINE(cert-err58-cpp)
static const std::array<std::string, 6> FILL_SAM_TAGS = {"AS", "XS", "XT", "XA", "SA", "MD"};

ReadCollector::ReadCollector(Params params) : mParams(std::move(params)) {
  using hts::Extractor;
  using hts::Alignment::Fields::AUX_RGAUX;

  mSampleList = MakeSampleList(mParams);
  const auto no_ctgcheck = mParams.mNoCtgCheck;

  for (const auto& sinfo : mSampleList) {
    auto extractor = std::make_unique<Extractor>(sinfo.Path(), mParams.mRefPath, AUX_RGAUX, FILL_SAM_TAGS, no_ctgcheck);
    mExtractors.emplace(sinfo, std::move(extractor));
  }
}

auto ReadCollector::CollectRegionResult(const Region& region) -> Result {
  std::vector<Read> sample_reads;
  absl::flat_hash_map<std::string, hts::Alignment::MateInfo> expected_mates;
  const auto max_sample_cov = mParams.mMaxWinCov / static_cast<f64>(mSampleList.size());

  for (auto& sinfo : mSampleList) {
    u64 num_reads = 0;
    u64 num_bases = 0;
    expected_mates.clear();

    const AlnAndRefPaths aln_refs{sinfo.Path(), mParams.mRefPath};
    const auto sample_name = std::string(sinfo.SampleName());
    const auto is_tumor_sample = sinfo.TagKind() == cbdg::Label::TUMOR;
    const auto est_sample_cov = EstimateCoverage(sinfo, region);

    const auto pct_to_sample = est_sample_cov > max_sample_cov ? ((max_sample_cov * 100.0) / est_sample_cov) : 100.0;
    mDownsampler.SetPercentToSample(pct_to_sample);

    auto& extractor = mExtractors.at(sinfo);
    extractor->SetRegionToExtract(region.ToSamtoolsRegion());

    for (const auto& aln : *extractor) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (FailsTier1Check(aln)) continue;
      if (is_tumor_sample && FailsTier2Check(aln)) continue;
      if (!mDownsampler.ShouldSample()) continue;
      // NOLINTEND(readability-braces-around-statements)

      num_reads += 1;
      num_bases += aln.Length();
      sample_reads.emplace_back(aln, sample_name, sinfo.TagKind());

      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!mParams.mExtractPairs) continue;

      // First check if we already saw both mates in the same window
      const auto mate_itr = expected_mates.find(aln.QnameView());
      if (mate_itr != expected_mates.end()) {
        expected_mates.erase(mate_itr);
        continue;
      }

      // NOLINTBEGIN(readability-braces-around-statements)
      if (aln.Flag().IsMateUnmapped() || FailsTier2Check(aln)) continue;
      if (aln.Flag().IsMappedProperPair() && !aln.HasTag("SA")) continue;
      // NOLINTEND(readability-braces-around-statements)

      auto [itr, newly_added] = expected_mates.try_emplace(aln.QnameView(), aln.MateLocation());
      // If not newly emplaced, then we already read both pairs
      // NOLINTNEXTLINE(readability-braces-around-statements)
      if (!newly_added) expected_mates.erase(itr);
    }

    if (!mParams.mExtractPairs || expected_mates.empty()) {
      sinfo.SetNumReads(num_reads);
      sinfo.SetNumBases(num_bases);
      sinfo.CalculateMeanCov(region.Length());
      continue;
    }

    auto rev_mate_regions = RevSortMateRegions(expected_mates);
    while (!expected_mates.empty()) {
      const auto& [rname, minfo] = rev_mate_regions.back();
      if (!expected_mates.contains(rname)) {
        rev_mate_regions.pop_back();
        continue;
      }

      const auto mate_reg_spec = MakeRegSpec(minfo, extractor.get());
      extractor->SetRegionToExtract(mate_reg_spec);

      for (const auto& aln : *extractor) {
        const auto itr = expected_mates.find(aln.QnameView());
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (itr == expected_mates.end()) continue;

        num_reads += 1;
        num_bases += aln.Length();
        sample_reads.emplace_back(aln, sample_name, sinfo.TagKind());
        expected_mates.erase(itr);
      }

      rev_mate_regions.pop_back();
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

auto ReadCollector::IsActiveRegion(const Params& params, const Region& region) -> bool {
  absl::btree_map<u32, u32> mismatches;  // genome position -> number of mismatches at position
  absl::btree_map<u32, u32> insertions;  // genome position -> number of insertions at position
  absl::btree_map<u32, u32> deletions;   // genome position -> number of deletions at position
  absl::btree_map<u32, u32> softclips;   // genome position -> number of softclips at position
  std::vector<u32> genome_positions;     // softclip genome positions for single alignment
  static const auto is_normal = [](const SampleInfo& sinfo) -> bool { return sinfo.TagKind() == cbdg::Label::NORMAL; };

  const auto sample_list = MakeSampleList(params);
  const auto is_germline_mode = std::ranges::all_of(sample_list, is_normal);

  for (const auto& sinfo : sample_list) {
    genome_positions.clear();
    mismatches.clear();
    insertions.clear();
    deletions.clear();
    softclips.clear();

    using hts::Alignment::Fields::AUX_RGAUX;
    const auto is_tumor_sample = sinfo.TagKind() == cbdg::Label::TUMOR;
    hts::Extractor extractor(sinfo.Path(), params.mRefPath, AUX_RGAUX, FILL_SAM_TAGS, params.mNoCtgCheck);
    extractor.SetRegionToExtract(region.ToSamtoolsRegion());

    for (const auto& aln : extractor) {
      // NOLINTBEGIN(readability-braces-around-statements)
      if (FailsTier1Check(aln) || aln.Flag().IsUnmapped()) continue;
      if ((is_tumor_sample || is_germline_mode) && FailsTier2Check(aln)) continue;
      // NOLINTEND(readability-braces-around-statements)

      if (aln.HasTag("MD")) {
        const auto md_tag = aln.GetTag<std::string_view>("MD");
        ParseMd(md_tag.value(), aln.QualView(), aln.StartPos0(), &mismatches);
      }

      const auto cigar_units = aln.CigarData();
      auto curr_genome_pos = static_cast<u32>(aln.StartPos0());
      bool has_soft_clip = false;

      // lambda function to increment the counter for `genome_positions`
      static const auto increment_genome_pos = [](CountMap& counts, u32 genome_pos) {
        auto [itr, newly_added] = counts.try_emplace(genome_pos, 1);
        // NOLINTNEXTLINE(readability-braces-around-statements)
        if (!newly_added) itr->second++;
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
        const auto handle_clips = [&softclips](const auto gpos) { increment_genome_pos(softclips, gpos); };
        std::ranges::for_each(genome_positions, handle_clips);
        LANCET_ASSERT(!softclips.empty())
      }

      // lambda function to check if `map` has evidence of mutation >= 2 reads
      static const auto count_gt2 = [](const CountMap& counts) {
        return std::any_of(counts.cbegin(), counts.cend(), [](const auto& pair) { return pair.second >= 2; });
      };

      // if evidence of SNV/insertion/deletion/softclip is found in window, mark region as active
      if (count_gt2(mismatches) || count_gt2(insertions) || count_gt2(deletions) || count_gt2(softclips)) {
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
  using hts::Alignment::Fields::AUX_RGAUX;
  const auto need_pairs = mParams.mExtractPairs;
  const auto need_tier2 = sinfo.TagKind() == cbdg::Label::TUMOR;

  hts::Extractor extractor(sinfo.Path(), mParams.mRefPath, AUX_RGAUX, FILL_SAM_TAGS, true);
  extractor.SetRegionToExtract(region.ToSamtoolsRegion());

  static const auto summer = [&need_tier2, &need_pairs, &region](const u64 sum, const hts::Alignment& aln) -> u64 {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (FailsTier1Check(aln) || (need_tier2 && FailsTier2Check(aln))) return sum;
    return (need_pairs && !aln.MateOverlapsRegion(region)) ? sum + aln.Length() + aln.Length() : sum + aln.Length();
  };

  const u64 num_bases = std::accumulate(extractor.begin(), extractor.end(), 0, summer);
  return static_cast<f64>(num_bases) / static_cast<f64>(region.Length());
}

auto ReadCollector::FailsTier1Check(const hts::Alignment& aln) -> bool {
  const auto bflag = aln.Flag();
  return bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsSecondary() || (bflag.IsMapped() && aln.MapQual() == 0);
}

auto ReadCollector::FailsTier2Check(const hts::Alignment& aln) -> bool {
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
  // XA -- BWA (Illumina): alternative hits; format: (chr,pos,CIGAR,NM;)
  return aln.HasTag("XT") || aln.HasTag("XA");
}

auto ReadCollector::MakeSampleList(const Params& params) -> std::vector<SampleInfo> {
  std::vector<SampleInfo> results;
  results.reserve(params.mNormalPaths.size() + params.mTumorPaths.size());
  const hts::Reference ref(params.mRefPath);

  // Add normal samples
  std::ranges::transform(params.mNormalPaths, std::back_inserter(results), [&ref](const std::filesystem::path& fpath) {
    const hts::Extractor extractor(fpath, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
    auto result = SampleInfo(extractor.SampleName(), fpath, cbdg::Label::NORMAL);
    return result;
  });

  // Add tumor samples
  std::ranges::transform(params.mTumorPaths, std::back_inserter(results), [&ref](const std::filesystem::path& fpath) {
    const hts::Extractor extractor(fpath, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
    auto result = SampleInfo(extractor.SampleName(), fpath, cbdg::Label::TUMOR);
    return result;
  });

  std::ranges::sort(results, std::less<SampleInfo>{});
  return results;
}

auto ReadCollector::RevSortMateRegions(const MateRegionsMap& data) -> std::vector<MateNameAndLocation> {
  std::vector<MateNameAndLocation> results(data.cbegin(), data.cend());
  std::ranges::sort(results, [](const MateNameAndLocation& lhs, const MateNameAndLocation& rhs) -> bool {
    return (lhs.second.mChromIndex != rhs.second.mChromIndex) ? lhs.second.mChromIndex > rhs.second.mChromIndex
                                                              : lhs.second.mMateStartPos0 > rhs.second.mMateStartPos0;
  });

  return results;
}

auto ReadCollector::MakeRegSpec(const hts::Alignment::MateInfo& info, const hts::Extractor* ext) -> std::string {
  const auto mate_chrom = ext->ChromName(info.mChromIndex);
  const auto mate_pos1 = info.mMateStartPos0 + 1;
  const auto colon_in_mate_chrom = mate_chrom.find(':') != std::string::npos;
  return colon_in_mate_chrom ? fmt::format("{{{}}}:{}-{}", mate_chrom, mate_pos1, mate_pos1)
                             : fmt::format("{}:{}-{}", mate_chrom, mate_pos1, mate_pos1);
}

}  // namespace lancet::core

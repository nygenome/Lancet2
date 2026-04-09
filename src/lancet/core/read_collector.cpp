#include "lancet/core/read_collector.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/cigar_unit.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/strings/ascii.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <random>
#include <spdlog/fmt/bundled/format.h>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdlib>

using CountMap = absl::flat_hash_map<u32, u32>;

namespace {

inline auto ParseMd(std::string_view md_val, absl::Span<u8 const> quals, i64 const start,
                    CountMap* result) -> bool {
  if (start < 0) return false;

  std::string token;
  token.reserve(md_val.length());
  auto genome_pos = static_cast<u32>(start);

  for (auto const& character : md_val) {
    if (absl::ascii_isdigit(static_cast<unsigned char>(character))) {
      token += character;
      continue;
    }

    auto const step = token.empty() ? 0 : std::strtol(token.c_str(), nullptr, 10);
    genome_pos += static_cast<u32>(step);
    token.clear();

    auto const base_pos = static_cast<usize>(genome_pos - start);
    static constexpr u8 MIN_BASE_QUAL = 20;
    if (quals.at(base_pos) < MIN_BASE_QUAL) continue;

    auto const base = absl::ascii_toupper(static_cast<unsigned char>(character));
    if (base == 'A' || base == 'C' || base == 'T' || base == 'G') {
      if (++(*result)[genome_pos] == 2) {
        return true;
      }
    }
  }

  return false;
}

}  // namespace

namespace lancet::core {

// ---------------------------------------------------------------------------
// Qname hashing utility
// ---------------------------------------------------------------------------

auto ReadCollector::HashQname(std::string_view qname) -> u64 {
  return absl::HashOf(qname);
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

ReadCollector::ReadCollector(Params params) : mParams(std::move(params)) {
  using hts::Extractor;
  using hts::Alignment::Fields::AUX_RGAUX;

  mSampleList = MakeSampleList(mParams);
  auto const no_ctgcheck = mParams.mNoCtgCheck;

  auto const sam_tags = mParams.mExtractPairs ? std::vector<std::string>{"SA", "AS", "XS"}
                                              : std::vector<std::string>{"AS", "XS"};
  static auto const IS_TUMOR = [](SampleInfo const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::TUMOR;
  };
  static auto const IS_NORMAL = [](SampleInfo const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::NORMAL;
  };
  mIsTumorNormalMode =
      std::ranges::any_of(mSampleList, IS_TUMOR) && std::ranges::any_of(mSampleList, IS_NORMAL);

  for (auto const& sinfo : mSampleList) {
    auto extractor = std::make_unique<Extractor>(sinfo.Path(), hts::Reference(mParams.mRefPath),
                                                 AUX_RGAUX, sam_tags, no_ctgcheck);
    mExtractors.emplace(sinfo, std::move(extractor));
  }
}

// ---------------------------------------------------------------------------
// CollectRegionResult: Two-pass uniform paired downsampling
//
// Pass 1 (Profile): Iterate the region with zero-copy alignment. Collect
//         per-read statistics (total/pass counts), build a set of unique
//         qname hashes, and track out-of-region mate candidates. Then
//         mathematically compute which template qnames to keep, ensuring
//         both mates of a pair are symmetrically accepted or rejected.
//
// Pass 2 (Extract): Re-iterate the region. Only reads whose qname hash
//         is in the keep_qnames set trigger deep extraction (BuildSequence,
//         BuildQualities via cbdg::Read construction).
//
// Pass 3 (Mates):   For kept reads whose mates fall outside the region,
//         fetch those mates using targeted single-position queries.
// ---------------------------------------------------------------------------

// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)  // TODO(lancet): refactor to reduce function size
auto ReadCollector::CollectRegionResult(Region const& region) -> Result {
  std::vector<Read> sampled_reads;
  auto const max_sample_bases = mParams.mMaxSampleCov * static_cast<f64>(region.Length());
  auto const region_spec = region.ToSamtoolsRegion();

  for (auto& sinfo : mSampleList) {
    auto& extractor = mExtractors.at(sinfo);
    auto const sample_name = std::string(sinfo.SampleName());

    // -----------------------------------------------------------------------
    // Pass 1: Profile & Downsample Math (zero-copy, no string allocations)
    // -----------------------------------------------------------------------
    u64 num_pass_reads = 0;
    u64 num_pass_bases = 0;

    // All unique qname hashes from passing reads (for downsampling)
    std::vector<u64> pass_qname_hashes;
    // Track out-of-region mate locations keyed by qname hash
    MateRegionsMap expected_mates;
    // Track qname hashes already seen in-region (for mate dedup)
    absl::flat_hash_set<u64> seen_in_region;

    extractor->SetRegionToExtract(region_spec);
    for (auto const& aln : *extractor) {
      auto const bflag = aln.Flag();
      if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() == 0)
        continue;

      auto const qhash = HashQname(aln.QnameView());
      bool const passes_filters = aln.MapQual() >= 20;
      if (passes_filters) {
        num_pass_reads += 1;
        num_pass_bases += aln.Length();
        pass_qname_hashes.push_back(qhash);
      }

      // Track out-of-region mates if pair extraction is enabled
      if (mParams.mExtractPairs) {
        // Check if we already saw the mate in-region
        if (seen_in_region.contains(qhash)) {
          expected_mates.erase(qhash);
        } else {
          seen_in_region.insert(qhash);

          if (!aln.Flag().IsMateUnmapped() &&
              (!aln.Flag().IsMappedProperPair() || aln.HasTag("SA"))) {
            expected_mates.try_emplace(qhash, aln.MateLocation());
          }
        }
      }
    }

    // Compute uniform paired downsampling threshold
    auto const bases_per_read =
        num_pass_reads > 0 ? static_cast<f64>(num_pass_bases) / static_cast<f64>(num_pass_reads)
                           : 1.0;
    auto const max_reads_to_sample = static_cast<u64>(std::ceil(max_sample_bases / bases_per_read));
    auto const sampled_read_count =
        num_pass_reads > max_reads_to_sample ? max_reads_to_sample : num_pass_reads;

    // Shuffle unique qname hashes and select the first `sampled_read_count` entries.
    // This guarantees that if Mate1 is accepted, Mate2 is symmetrically accepted.
    std::shuffle(pass_qname_hashes.begin(), pass_qname_hashes.end(),
                 std::default_random_engine(0));  // NOLINT(cert-msc51-cpp)
    absl::flat_hash_set<u64> const keep_qnames(pass_qname_hashes.begin(),
                                               pass_qname_hashes.begin() +
                                                   static_cast<i64>(sampled_read_count));

    // -----------------------------------------------------------------------
    // Pass 2: Deep Copy & Object Emplacement (only for kept reads)
    // -----------------------------------------------------------------------
    u64 sampled_base_count = 0;
    extractor->SetRegionToExtract(region_spec);
    for (auto const& aln : *extractor) {
      auto const bflag = aln.Flag();
      if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() == 0)
        continue;
      if (!keep_qnames.contains(HashQname(aln.QnameView()))) continue;

      // Only kept reads trigger BuildSequence/BuildQualities via the Read constructor
      sampled_reads.emplace_back(aln, sample_name, sinfo.TagKind());
      sampled_base_count += sampled_reads.back().Length();
    }

    // -----------------------------------------------------------------------
    // Pass 3: Recapture out-of-region mates for kept reads
    // -----------------------------------------------------------------------
    if (!expected_mates.empty() && mParams.mExtractPairs) {
      // Remove mates whose qnames were not kept during downsampling
      for (auto it = expected_mates.begin(); it != expected_mates.end();) {
        if (!keep_qnames.contains(it->first)) {
          expected_mates.erase(it++);
        } else {
          ++it;
        }
      }

      auto rev_mate_regions = RevSortMateRegions(expected_mates);
      while (!expected_mates.empty()) {
        auto const& [qhash, minfo] = rev_mate_regions.back();
        if (!expected_mates.contains(qhash)) {
          rev_mate_regions.pop_back();
          continue;
        }

        auto const mate_reg_spec = MakeRegSpec(minfo, extractor.get());
        extractor->SetRegionToExtract(mate_reg_spec);

        for (auto const& aln : *extractor) {
          auto const mate_qhash = HashQname(aln.QnameView());
          auto const itr = expected_mates.find(mate_qhash);
          if (itr == expected_mates.end()) continue;

          sampled_reads.emplace_back(aln, sample_name, sinfo.TagKind());
          sampled_base_count += sampled_reads.back().Length();
          expected_mates.erase(itr);
        }

        rev_mate_regions.pop_back();
      }
    }

    // Update sample statistics
    sinfo.SetNumSampledReads(sampled_read_count);
    sinfo.SetNumSampledBases(sampled_base_count);
  }

  std::ranges::sort(sampled_reads, [](Read const& lhs, Read const& rhs) -> bool {
    if (lhs.PassesAlnFilters() != rhs.PassesAlnFilters()) {
      return static_cast<int>(lhs.PassesAlnFilters()) > static_cast<int>(rhs.PassesAlnFilters());
    }
    if (lhs.TagKind() != rhs.TagKind())
      return static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind());
    if (lhs.SampleName() != rhs.SampleName()) return lhs.SampleName() < rhs.SampleName();
    if (lhs.QnameView() != rhs.QnameView()) return lhs.QnameView() < rhs.QnameView();
    if (lhs.ChromIndex() != rhs.ChromIndex()) return lhs.ChromIndex() < rhs.ChromIndex();
    return lhs.StartPos0() < rhs.StartPos0();
  });

  return {.mSampleReads = std::move(sampled_reads), .mSampleList = mSampleList};
}

// ---------------------------------------------------------------------------
// IsActiveRegion: uses zero-copy alignment for rapid region evaluation.
// BuildQualities() is called on-demand only for reads that have MD tags.
// ---------------------------------------------------------------------------

// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)  // TODO(lancet): refactor to reduce function size
auto ReadCollector::IsActiveRegion(Params const& params, Region const& region) -> bool {
  CountMap mismatches;                // genome position -> number of mismatches at position
  CountMap insertions;                // genome position -> number of insertions at position
  CountMap deletions;                 // genome position -> number of deletions at position
  CountMap softclips;                 // genome position -> number of softclips at position
  std::vector<u32> genome_positions;  // softclip genome positions for single alignment

  auto const sample_list = MakeSampleList(params);
  for (auto const& sinfo : sample_list) {
    genome_positions.clear();
    mismatches.clear();
    insertions.clear();
    deletions.clear();
    softclips.clear();

    using hts::Alignment::Fields::AUX_RGAUX;
    hts::Extractor extractor(sinfo.Path(), hts::Reference(params.mRefPath), AUX_RGAUX,
                             {"MD", "AS", "XS"}, params.mNoCtgCheck);
    extractor.SetRegionToExtract(region.ToSamtoolsRegion());

    for (auto const& aln : extractor) {
      auto const bflag = aln.Flag();
      if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() == 0)
        continue;

      if (aln.HasTag("MD")) {
        auto const md_tag = aln.GetTag<std::string_view>("MD");
        // BuildQualities performs on-demand deep copy of quality values
        auto const quals = aln.BuildQualities();
        if (ParseMd(md_tag.value(), absl::MakeConstSpan(quals), aln.StartPos0(), &mismatches)) {
          return true;
        }
      }

      auto const cigar_units = aln.CigarData();
      auto curr_genome_pos = static_cast<u32>(aln.StartPos0());
      bool has_soft_clip = false;

      // lambda function to increment the counter for `genome_positions`
      static auto const INCREMENT_GENOME_POS = [](CountMap& counts, u32 genome_pos) -> bool {
        return ++counts[genome_pos] == 2;
      };

      for (auto const& cig_unit : cigar_units) {
        if (cig_unit.ConsumesReference()) {
          curr_genome_pos += cig_unit.Length();
        }

        switch (cig_unit.Operation()) {
          case hts::CigarOp::INSERTION:
            if (INCREMENT_GENOME_POS(insertions, curr_genome_pos)) {
              return true;
            }
            break;
          case hts::CigarOp::DELETION:
            if (INCREMENT_GENOME_POS(deletions, curr_genome_pos)) {
              return true;
            }
            break;
          case hts::CigarOp::SEQUENCE_MISMATCH:
            if (INCREMENT_GENOME_POS(mismatches, curr_genome_pos)) {
              return true;
            }
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
        for (auto const gpos : genome_positions) {
          if (INCREMENT_GENOME_POS(softclips, gpos)) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

// ---------------------------------------------------------------------------
// Utility methods
// ---------------------------------------------------------------------------

auto ReadCollector::BuildSampleNameList(Params const& params) -> std::vector<std::string> {
  auto const sinfo_list = MakeSampleList(params);
  std::vector<std::string> results;
  results.reserve(sinfo_list.size());
  std::ranges::transform(
      sinfo_list, std::back_inserter(results),
      [](SampleInfo const& item) -> std::string { return std::string(item.SampleName()); });
  return results;
}

auto ReadCollector::MakeSampleList(Params const& params) -> std::vector<SampleInfo> {
  std::vector<SampleInfo> results;
  results.reserve(params.mNormalPaths.size() + params.mTumorPaths.size());
  hts::Reference const ref(params.mRefPath);

  // Add normal samples
  std::ranges::transform(
      params.mNormalPaths, std::back_inserter(results),
      [&ref](std::filesystem::path const& fpath) -> SampleInfo {
        hts::Extractor const extractor(fpath, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
        auto result = SampleInfo(extractor.SampleName(), fpath, cbdg::Label::NORMAL);
        return result;
      });

  // Add tumor samples
  std::ranges::transform(
      params.mTumorPaths, std::back_inserter(results),
      [&ref](std::filesystem::path const& fpath) -> SampleInfo {
        hts::Extractor const extractor(fpath, ref, hts::Alignment::Fields::CORE_QNAME, {}, true);
        auto result = SampleInfo(extractor.SampleName(), fpath, cbdg::Label::TUMOR);
        return result;
      });

  std::ranges::sort(results, std::less<SampleInfo>{});
  return results;
}

auto ReadCollector::RevSortMateRegions(MateRegionsMap const& data)
    -> std::vector<MateHashAndLocation> {
  std::vector<MateHashAndLocation> results(data.cbegin(), data.cend());
  std::ranges::sort(results,
                    [](MateHashAndLocation const& lhs, MateHashAndLocation const& rhs) -> bool {
                      return (lhs.second.mChromIndex != rhs.second.mChromIndex)
                                 ? lhs.second.mChromIndex > rhs.second.mChromIndex
                                 : lhs.second.mMateStartPos0 > rhs.second.mMateStartPos0;
                    });

  return results;
}

auto ReadCollector::MakeRegSpec(hts::Alignment::MateInfo const& info, hts::Extractor const* ext)
    -> std::string {
  auto const mate_chrom = ext->ChromName(info.mChromIndex);
  auto const mate_pos1 = info.mMateStartPos0 + 1;
  auto const colon_in_mate_chrom = mate_chrom.find(':') != std::string::npos;
  return colon_in_mate_chrom ? fmt::format("{{{}}}:{}-{}", mate_chrom, mate_pos1, mate_pos1)
                             : fmt::format("{}:{}-{}", mate_chrom, mate_pos1, mate_pos1);
}

}  // namespace lancet::core

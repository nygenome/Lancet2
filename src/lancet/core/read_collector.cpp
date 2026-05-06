#include "lancet/core/read_collector.h"

#include "lancet/base/types.h"
#include "lancet/cbdg/label.h"
#include "lancet/cbdg/read.h"
#include "lancet/core/sample_info.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/extractor.h"
#include "lancet/hts/mate_info.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"

#include <algorithm>
#include <array>
#include <functional>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdlib>

namespace {

// ============================================================================
// CompareReadsByPriority — 6-key composite comparator for deterministic
// read ordering. Ensures identical output regardless of HTSlib iteration
// order across runs.
//
// Priority: filter-pass status > sample tag > sample name > qname > chrom > position.
// ============================================================================
auto CompareReadsByPriority(lancet::cbdg::Read const& lhs, lancet::cbdg::Read const& rhs) -> bool {
  if (lhs.PassesAlnFilters() != rhs.PassesAlnFilters()) {
    return static_cast<int>(lhs.PassesAlnFilters()) > static_cast<int>(rhs.PassesAlnFilters());
  }
  if (lhs.TagKind() != rhs.TagKind()) {
    return static_cast<u8>(lhs.TagKind()) < static_cast<u8>(rhs.TagKind());
  }
  if (lhs.SampleName() != rhs.SampleName()) return lhs.SampleName() < rhs.SampleName();
  if (lhs.QnameView() != rhs.QnameView()) return lhs.QnameView() < rhs.QnameView();
  if (lhs.ChromIndex() != rhs.ChromIndex()) return lhs.ChromIndex() < rhs.ChromIndex();
  return lhs.StartPos0() < rhs.StartPos0();
}

}  // namespace

namespace lancet::core {

// ============================================================================
// Qname hashing utility
// ============================================================================
auto ReadCollector::HashQname(std::string_view qname) -> u64 {
  return absl::HashOf(qname);
}

// ============================================================================
// Constructor
// ============================================================================
ReadCollector::ReadCollector(Params params, absl::Span<SampleInfo const> sample_list)
    : mParams(std::move(params)), mSampleList(sample_list.begin(), sample_list.end()) {
  using hts::Extractor;
  using hts::Alignment::Fields::AUX_RGAUX;

  auto const no_ctgcheck = mParams.mNoCtgCheck;

  // Always request MD (for IsActiveRegion prescan), conditionally SA (for mate pairs).
  // HTSlib lazily parses AUX fields — adding MD has negligible per-read cost.
  static std::array<std::string, 1> const MD_TAG{"MD"};
  static std::array<std::string, 2> const MD_SA_TAGS{"MD", "SA"};
  auto const sam_tags =
      mParams.mExtractPairs ? absl::MakeConstSpan(MD_SA_TAGS) : absl::MakeConstSpan(MD_TAG);

  static auto const IS_CASE = [](SampleInfo const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::CASE;
  };
  static auto const IS_CTRL = [](SampleInfo const& sinfo) -> bool {
    return sinfo.TagKind() == cbdg::Label::CTRL;
  };
  mIsCaseCtrlMode =
      std::ranges::any_of(mSampleList, IS_CASE) && std::ranges::any_of(mSampleList, IS_CTRL);

  for (auto const& sinfo : mSampleList) {
    auto extractor = std::make_unique<Extractor>(sinfo.Path(), hts::Reference(mParams.mRefPath),
                                                 AUX_RGAUX, sam_tags, no_ctgcheck);
    mExtractors.emplace(sinfo, std::move(extractor));
  }
}

// ============================================================================
// CollectRegionResult: orchestrator for three-pass paired downsampling.
//
// Pass 1 (Profile): zero-copy profiling + deterministic downsampling.
// Pass 2 (Extract): deep-copy only kept reads into mSampledReads.
// Pass 3 (Mates):   fetch out-of-region mates for kept reads.
// ============================================================================
auto ReadCollector::CollectRegionResult(Region const& region) -> Result {
  mSampledReads.clear();
  auto const max_sample_bases = mParams.mMaxSampleCov * static_cast<f64>(region.Length());
  auto const region_spec = region.ToSamtoolsRegion();

  for (auto& sinfo : mSampleList) {
    auto& extractor = mExtractors.at(sinfo);
    mSampledBaseCount = 0;

    auto profile = ProfileAndDownsample(*extractor, region_spec, max_sample_bases);
    ExtractKeptReads(*extractor, region_spec, profile.mKeepQnames, sinfo);

    if (!profile.mExpectedMates.empty() && mParams.mExtractPairs) {
      RecaptureMates(*extractor, profile.mKeepQnames, profile.mExpectedMates, sinfo);
    }

    sinfo.SetNumSampledReads(profile.mSampledReadCount);
    sinfo.SetNumSampledBases(mSampledBaseCount);
  }

  std::ranges::sort(mSampledReads, CompareReadsByPriority);
  return {.mSampleReads = std::move(mSampledReads), .mSampleList = mSampleList};
}

// ============================================================================
// Pass 1: Profile & Downsample Math (zero-copy, no string allocations)
//
// Iterates all alignments in the region for a single sample. Counts
// passing reads/bases, builds a qname hash set, and tracks out-of-region
// mate locations. Then shuffles the qname hashes and keeps only enough
// to satisfy the coverage cap. Both mates of a pair are symmetrically
// accepted or rejected because downsampling operates on qname hashes.
// ============================================================================
auto ReadCollector::ProfileAndDownsample(hts::Extractor& extractor, std::string const& region_spec,
                                         f64 const max_sample_bases) const -> ProfileResult {
  u64 num_pass_reads = 0;
  u64 num_pass_bases = 0;

  std::vector<u64> pass_qname_hashes;
  MateRegionsMap expected_mates;
  absl::flat_hash_set<u64> seen_in_region;

  extractor.SetRegionToExtract(region_spec);
  for (auto const& aln : extractor) {
    auto const bflag = aln.Flag();
    if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() < 20) {
      continue;
    }

    auto const qhash = HashQname(aln.QnameView());
    num_pass_reads += 1;
    num_pass_bases += aln.Length();
    pass_qname_hashes.push_back(qhash);

    if (!mParams.mExtractPairs) continue;

    // Both mates already seen in-region — no out-of-region fetch needed
    if (seen_in_region.contains(qhash)) {
      expected_mates.erase(qhash);
      continue;
    }

    seen_in_region.insert(qhash);
    // Track mates with discordant mapping or supplementary alignments (SA tag)
    auto const is_mate_mapped = aln.Flag().IsMateMapped();
    auto const is_discordant = !aln.Flag().IsMappedProperPair();
    auto const has_supplementary = aln.HasTag("SA");
    if (is_mate_mapped && (is_discordant || has_supplementary)) {
      expected_mates.try_emplace(qhash, aln.MateLocation());
    }
  }

  // Compute the maximum number of reads to keep for coverage-capped downsampling.
  // Mean read length converts the base-budget (max_sample_bases) into a read count.
  // std::max guards against zero-read regions to avoid division by zero (branchless).
  // std::min caps at observed count so we never "upsample" beyond what exists.
  auto const denominator = static_cast<f64>(std::max(num_pass_reads, u64{1}));
  auto const bases_per_read = static_cast<f64>(num_pass_bases) / denominator;
  auto const max_reads_to_sample = static_cast<u64>(std::ceil(max_sample_bases / bases_per_read));
  auto const sampled_read_count = std::min(num_pass_reads, max_reads_to_sample);

  // Shuffle and select the first `sampled_read_count` entries.
  // Fixed seed=0 is intentional: variant calls MUST be deterministic across
  // runs on identical inputs. A random seed would cause non-reproducible VCF
  // output, violating the pipeline's reproducibility guarantee.
  // The engine is the C++-standard seed-stable Mersenne Twister; its output
  // sequence is fixed by the C++ standard so the same seed produces the
  // same shuffle order on any conforming compiler.
  // NOLINTBEGIN(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::shuffle(pass_qname_hashes.begin(), pass_qname_hashes.end(), std::mt19937_64{0});
  // NOLINTEND(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)

  auto const keep_qnames_end = pass_qname_hashes.begin() + static_cast<i64>(sampled_read_count);
  absl::flat_hash_set<u64> keep_qnames(pass_qname_hashes.begin(), keep_qnames_end);

  return {.mKeepQnames = std::move(keep_qnames),
          .mExpectedMates = std::move(expected_mates),
          .mSampledReadCount = sampled_read_count};
}

// ============================================================================
// Pass 2: Deep Copy & Object Emplacement (only for kept reads)
//
// Re-iterates the region. Only reads whose qname hash is in keep_qnames
// trigger deep extraction (BuildSequence, BuildQualities via Read ctor).
// ============================================================================
void ReadCollector::ExtractKeptReads(hts::Extractor& extractor, std::string const& region_spec,
                                     absl::flat_hash_set<u64> const& keep_qnames,
                                     SampleInfo const& sinfo) {
  auto const sample_name = std::string(sinfo.SampleName());
  extractor.SetRegionToExtract(region_spec);
  for (auto const& aln : extractor) {
    auto const bflag = aln.Flag();
    if (bflag.IsQcFail() || bflag.IsDuplicate() || bflag.IsUnmapped() || aln.MapQual() < 20) {
      continue;
    }

    if (!keep_qnames.contains(HashQname(aln.QnameView()))) continue;

    mSampledReads.emplace_back(aln, sample_name, sinfo.TagKind(), sinfo.SampleIndex());
    mSampledBaseCount += mSampledReads.back().Length();
  }
}

// ============================================================================
// Pass 3: Recapture out-of-region mates for kept reads
//
// Walks mate locations in reverse-sorted genomic order for sequential
// disk access. Removes mates whose qnames were not kept during downsampling
// before querying, to avoid unnecessary HTSlib seeks.
// ============================================================================
void ReadCollector::RecaptureMates(hts::Extractor& extractor,
                                   absl::flat_hash_set<u64> const& keep_qnames,
                                   MateRegionsMap& expected_mates, SampleInfo const& sinfo) {
  // Remove mates whose qnames are not in the keep set (rejected by downsampling
  // or already captured during the in-region linear scan)
  absl::erase_if(expected_mates, [&keep_qnames](auto const& entry) -> bool {
    return !keep_qnames.contains(entry.first);
  });

  auto const sample_name = std::string(sinfo.SampleName());
  // Reverse-sorted so pop_back() yields ascending genomic order (sequential BAM I/O)
  auto rev_mate_regions = ReverseSortMateRegions(expected_mates);
  while (!expected_mates.empty()) {
    auto const& [qhash, minfo] = rev_mate_regions.back();
    if (!expected_mates.contains(qhash)) {
      rev_mate_regions.pop_back();
      continue;
    }

    auto const mate_reg_spec = MakeRegionSpec(minfo, &extractor);
    extractor.SetRegionToExtract(mate_reg_spec);

    for (auto const& aln : extractor) {
      auto const mate_qhash = HashQname(aln.QnameView());
      auto const itr = expected_mates.find(mate_qhash);
      if (itr == expected_mates.end()) continue;

      mSampledReads.emplace_back(aln, sample_name, sinfo.TagKind(), sinfo.SampleIndex());
      mSampledBaseCount += mSampledReads.back().Length();
      expected_mates.erase(itr);
    }

    rev_mate_regions.pop_back();
  }
}

// ============================================================================
// ReverseSortMateRegions — descending-order sort for stack-based consumption
//
// Returns mate locations sorted by descending chrom index, then descending
// start position. RecaptureMates consumes via pop_back(), which yields
// ascending genomic order — matching BAM/CRAM coordinate sort for sequential
// disk I/O. Using pop_back() is O(1), unlike front-erasure which is O(n).
// ============================================================================
auto ReadCollector::ReverseSortMateRegions(MateRegionsMap const& data)
    -> std::vector<MateHashAndLocation> {
  // Descending comparator: higher chrom index first,
  // then higher start position within the same chrom.
  static auto const COMPARE_DESCENDING = [](MateHashAndLocation const& lhs,
                                            MateHashAndLocation const& rhs) -> bool {
    if (lhs.second.mChromIndex != rhs.second.mChromIndex) {
      return lhs.second.mChromIndex > rhs.second.mChromIndex;
    }
    return lhs.second.mMateStartPos0 > rhs.second.mMateStartPos0;
  };

  std::vector<MateHashAndLocation> results(data.cbegin(), data.cend());
  std::ranges::sort(results, COMPARE_DESCENDING);
  return results;
}

auto ReadCollector::MakeRegionSpec(hts::MateInfo const& info, hts::Extractor const* ext)
    -> std::string {
  auto const mate_chrom = ext->ChromName(info.mChromIndex);
  auto const mate_pos1 = info.mMateStartPos0 + 1;
  // HTSlib region syntax: chroms with ':' in their name (e.g. HLA contigs)
  // must be brace-wrapped as {chrom}:pos-pos to avoid ambiguous parsing.
  auto const colon_in_mate_chrom = mate_chrom.find(':') != std::string::npos;
  return colon_in_mate_chrom ? fmt::format("{{{}}}:{}-{}", mate_chrom, mate_pos1, mate_pos1)
                             : fmt::format("{}:{}-{}", mate_chrom, mate_pos1, mate_pos1);
}

}  // namespace lancet::core

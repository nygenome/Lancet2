#include "lancet/hts/reference.h"

#include <algorithm>
#include <filesystem>
#include <optional>
#include <spdlog/fmt/bundled/format.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cstdlib>

extern "C" {
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/hts_log.h"
}

#include "lancet/base/types.h"

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "spdlog/fmt/bundled/core.h"

namespace lancet::hts {

Reference::Reference(std::filesystem::path reference) : mFastaPath(std::move(reference)) {
  hts_set_log_level(HTS_LOG_ERROR);
  mFastaIndex = FastaIndex(fai_load3(mFastaPath.c_str(), nullptr, nullptr, 0));
  if (mFastaIndex == nullptr) {
    auto const fname = mFastaPath.filename().string();
    auto const msg = fmt::format("Could not load mIdx for reference: {}", fname);
    throw std::runtime_error(msg);
  }

  // Set cache size enough to load all blocks from a chrom into memory i.e. ~268 MB
  static constexpr int DEFAULT_FAI_CACHE_SIZE = 1 << 28;
  fai_set_cache_size(mFastaIndex.get(), DEFAULT_FAI_CACHE_SIZE);

  auto const num_chroms = faidx_nseq(mFastaIndex.get());
  if (num_chroms <= 0) {
    auto const fname = mFastaPath.filename().string();
    auto const msg = fmt::format("No chromosomes found in reference: {}", fname);
    throw std::runtime_error(msg);
  }

  mChroms.reserve(static_cast<usize>(num_chroms));
  for (auto idx = 0; idx < num_chroms; ++idx) {
    auto const chrom_name = std::string(faidx_iseq(mFastaIndex.get(), idx));
    auto const chrom_length = faidx_seq_len64(mFastaIndex.get(), chrom_name.c_str());
    mChroms.emplace_back(Chrom(idx, chrom_name, chrom_length));
  }
}

auto Reference::ListChroms() const noexcept -> std::vector<Chrom> {
  return mChroms;
}

auto Reference::FindChromByName(std::string_view chrom_name) const noexcept
    -> absl::StatusOr<Chrom> {
  auto const itr = std::ranges::find_if(
      mChroms, [&chrom_name](Chrom const& chrom) -> bool { return chrom.Name() == chrom_name; });

  if (itr == mChroms.cend()) {
    auto const msg =
        fmt::format("Chrom {} not found in reference: {}", chrom_name, mFastaPath.string());
    return absl::Status(absl::StatusCode::kNotFound, msg);
  }

  return *itr;
}

auto Reference::FindChromByIndex(i64 chrom_index) const noexcept -> absl::StatusOr<Chrom> {
  auto const nchroms = mChroms.size();
  if (chrom_index < 0 || std::cmp_greater_equal(chrom_index, nchroms)) {
    return absl::Status(absl::StatusCode::kOutOfRange,
                        fmt::format("Index {} is out of range for reference with {} chromosomes",
                                    chrom_index, nchroms));
  }

  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  return mChroms[chrom_index];
}

auto Reference::MakeRegion(std::string const& chrom_name,
                           OneBasedClosedOptional const& interval) const -> Region {
  auto const matching_chrom = FindChromByName(chrom_name);
  if (!matching_chrom.ok()) {
    auto const msg =
        fmt::format("Chrom {} not found in reference: {}", chrom_name, mFastaPath.string());
    throw std::invalid_argument(msg);
  }

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  u64 const given_start = interval[0].value_or(1);
  u64 const given_end = interval[1].value_or(matching_chrom->Length());
  // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)

  if (given_start == 0 || given_end == 0) {
    throw std::invalid_argument("Expected 1-based co-ordinates for start and end positions");
  }

  if (given_start > matching_chrom->Length() || given_end > matching_chrom->Length()) {
    throw std::out_of_range("Expected start and end positions to be <= chromosome mVarLength");
  }

  if (given_start > given_end) {
    throw std::invalid_argument("Expected start position to be <= end position");
  }

  auto region_seq = FetchSeq(chrom_name, {given_start, given_end});
  return {
      matching_chrom->Index(), {given_start, given_end}, chrom_name.c_str(), std::move(region_seq)};
}

auto Reference::MakeRegion(ParseRegionResult const& parse_result) const -> Region {
  return MakeRegion(parse_result.mChromName, parse_result.mRegionSpan);
}

auto Reference::MakeRegion(char const* region_spec) const -> Region {
  auto const result = ParseRegion(region_spec);
  return MakeRegion(result.mChromName, result.mRegionSpan);
}

auto Reference::ParseRegion(char const* region_spec) const -> ParseRegionResult {
  int tid = 0;
  hts_pos_t begin_zero = 0;
  hts_pos_t end_one = 0;
  auto const* res = fai_parse_region(mFastaIndex.get(), region_spec, &tid, &begin_zero, &end_one,
                                     HTS_PARSE_ONE_COORD);
  if (res == nullptr || tid < 0 || begin_zero < 0 || end_one < (begin_zero + 1)) {
    auto const err_msg =
        fmt::format("Could not parse string as a samtools region: {}", region_spec);
    throw std::runtime_error(err_msg);
  }

  auto const adjust_result = fai_adjust_region(mFastaIndex.get(), tid, &begin_zero, &end_one);
  if (adjust_result < 0) {
    auto const err_msg =
        fmt::format("Could not adjust co-ordinates of samtools region: {}", region_spec);
    throw std::runtime_error(err_msg);
  }

  return ParseRegionResult{.mChromName = faidx_iseq(mFastaIndex.get(), tid),
                           .mRegionSpan = {begin_zero + 1, end_one}};
}

auto Reference::FetchSeq(std::string const& chrom, OneBasedClosedInterval const& full_intvl) const
    -> std::string {
  hts_pos_t parsed_len = 0;
  auto const [start_pos1, end_pos1] = full_intvl;
  char* raw_seq =
      faidx_fetch_seq64(mFastaIndex.get(), chrom.c_str(), static_cast<hts_pos_t>(start_pos1 - 1),
                        static_cast<hts_pos_t>(end_pos1 - 1), &parsed_len);
  if (raw_seq == nullptr) {
    auto const colon_in_chrom = chrom.find(':') != std::string::npos;
    auto const regspec = colon_in_chrom ? fmt::format("{{{}}}:{}-{}", chrom, start_pos1, end_pos1)
                                        : fmt::format("{}:{}-{}", chrom, start_pos1, end_pos1);
    auto const err_msg = fmt::format("Got zero mVarLength mDfltSeq from region: {}", regspec);
    throw std::runtime_error(err_msg);
  }

  auto const expected_length = end_pos1 - start_pos1 + 1;
  if (std::cmp_not_equal(parsed_len, expected_length)) {
    constexpr auto FORMAT_MSG = "Expected to get {} bases from region {}. Got {} bases instead";
    auto const colon_in_chrom = chrom.find(':') != std::string::npos;
    auto const regspec = colon_in_chrom ? fmt::format("{{{}}}:{}-{}", chrom, start_pos1, end_pos1)
                                        : fmt::format("{}:{}-{}", chrom, start_pos1, end_pos1);
    auto const err_msg = fmt::format(FORMAT_MSG, expected_length, regspec, parsed_len);
    throw std::runtime_error(err_msg);
  }

  auto const sequence_length = static_cast<usize>(parsed_len);
  std::string result_seq(sequence_length, 'N');
  for (usize base_idx = 0; base_idx < sequence_length; ++base_idx) {
    auto const base = absl::ascii_toupper(static_cast<unsigned char>(raw_seq[base_idx]));  // NOLINT
    switch (base) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        result_seq[base_idx] = base;
        break;
      default:
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        result_seq[base_idx] = 'N';
    }
  }

  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
  std::free(raw_seq);
  return std::move(result_seq);
}

auto Reference::Region::ToSamtoolsRegion() const -> std::string {
  auto const name_has_colon = mName.find(':') != std::string::npos;
  return name_has_colon ? fmt::format("{{{}}}:{}-{}", mName, mStart1, mEnd1)
                        : fmt::format("{}:{}-{}", mName, mStart1, mEnd1);
}

auto Reference::ParseRegionResult::Length() const noexcept -> usize {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
  if (mRegionSpan[0] == std::nullopt || mRegionSpan[1] == std::nullopt) {
    return 0;
  }

  return mRegionSpan[1].value_or(1) - mRegionSpan[0].value_or(0) + 1;
  // NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
}

}  // namespace lancet::hts

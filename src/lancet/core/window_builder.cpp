#include "lancet/core/window_builder.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "absl/types/span.h"
#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/core/window.h"
#include "spdlog/fmt/bundled/core.h"

namespace lancet::core {

WindowBuilder::WindowBuilder(const std::filesystem::path &ref_path, const Params &params)
    : mParams(params), mRefPtr(std::make_unique<hts::Reference>(ref_path)) {
  static constexpr usize DEFAULT_NUM_REGIONS_TO_ALLOCATE = 1024;
  mInputRegions.reserve(DEFAULT_NUM_REGIONS_TO_ALLOCATE);
}

void WindowBuilder::AddAllReferenceRegions() {
  static const auto should_exclude_chrom = [](std::string_view chrom) -> bool {
    return chrom == "MT" || chrom == "chrM" || absl::StartsWith(chrom, "GL") || absl::StartsWith(chrom, "chrUn") ||
           absl::StartsWith(chrom, "chrEBV") || absl::StartsWith(chrom, "HLA-") || absl::EndsWith(chrom, "_random") ||
           absl::EndsWith(chrom, "_alt") || absl::EndsWith(chrom, "_decoy");
  };

  const auto ref_chroms = mRefPtr->ListChroms();
  mInputRegions.reserve(mInputRegions.size() + ref_chroms.size());

  std::ranges::for_each(ref_chroms, [this](const hts::Reference::Chrom &chrom) {
    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (should_exclude_chrom(chrom.Name())) return;
    this->mInputRegions.emplace(ParseRegionResult{.mChromName = chrom.Name(), .mRegionSpan = {1, chrom.Length()}});
  });
}

void WindowBuilder::AddRegion(const std::string &region_spec) {
  mInputRegions.emplace(mRefPtr->ParseRegion(region_spec.c_str()));
}

void WindowBuilder::AddBatchRegions(absl::Span<const std::string> region_specs) {
  mInputRegions.reserve(mInputRegions.size() + region_specs.size());
  std::ranges::for_each(region_specs, [this](const std::string &spec) {
    this->mInputRegions.emplace(this->mRefPtr->ParseRegion(spec.c_str()));
  });
}

void WindowBuilder::AddBatchRegions(const std::filesystem::path &bed_file) {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (bed_file.empty() || !std::filesystem::exists(bed_file)) return;

  std::ifstream fhandle(bed_file, std::ios_base::in);
  std::string contents;

  if (fhandle) {
    fhandle.seekg(0, std::ifstream::end);
    const std::int64_t length = fhandle.tellg();
    fhandle.seekg(0, std::ifstream::beg);
    contents.resize(static_cast<std::size_t>(length), '\0');
    fhandle.read(contents.data(), length);
    fhandle.close();
  }

  usize line_num = 0;
  i64 region_start = 0;
  i64 region_end = 0;

  std::string curr_chrom;
  std::vector<std::string_view> tokens;
  tokens.reserve(3);

  for (const auto &line : absl::StrSplit(contents, absl::ByChar('\n'))) {
    line_num++;

    // NOLINTNEXTLINE(readability-braces-around-statements)
    if (line.starts_with('#') || line.empty()) continue;

    tokens = absl::StrSplit(line, absl::ByChar('\t'));

    if (tokens.size() != 3) {
      const auto msg = fmt::format("Invalid bed line with {} columns at line number {}", tokens.size(), line_num);
      throw std::runtime_error(msg);
    }

    if (!absl::SimpleAtoi(tokens[1], &region_start) || !absl::SimpleAtoi(tokens[2], &region_end)) {
      const auto msg = fmt::format("Could not parse line {} in bed: {}", line_num, bed_file.filename().string());
      throw std::runtime_error(msg);
    }

    curr_chrom = tokens[0];
    if (!mRefPtr->FindChromByName(curr_chrom).ok()) {
      const auto msg = fmt::format("Could not find chrom {} from bed file line {} in reference", tokens[0], line_num);
      throw std::runtime_error(msg);
    }

    mInputRegions.emplace(ParseRegionResult{.mChromName = curr_chrom, .mRegionSpan = {region_start, region_end}});
  }
}

auto WindowBuilder::StepSize(const Params &params) -> i64 {
  const auto val = (static_cast<f64>(100 - params.mPercentOverlap) / 100.0) * static_cast<f64>(params.mWindowLength);
  // round to ensure that steps always move in multiples of 100
  return static_cast<i64>(std::ceil(val / 100.0) * 100.0);
}

auto WindowBuilder::BuildWindows() const -> std::vector<WindowPtr> {
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (mInputRegions.empty()) return {};

  std::string rspec;
  const auto nregs = mInputRegions.size();
  const auto window_len = static_cast<i64>(mParams.mWindowLength);
  const auto pct_olap = static_cast<i64>(mParams.mPercentOverlap);
  LOG_INFO("Using {} input region(s) to build {}bp moving windows with {}% overlap", nregs, window_len, pct_olap)

  const auto step_size = StepSize(mParams);
  absl::flat_hash_set<WindowPtr> uniq_windows;

  for (ParseRegionResult region : mInputRegions) {
    PadInputRegion(region);
    const auto chrom = mRefPtr->FindChromByName(region.mChromName).value();

    if (region.Length() <= window_len) {
      auto wptr = std::make_shared<Window>(std::move(region), chrom, mRefPtr->FastaPath());
      uniq_windows.emplace(std::move(wptr));
      continue;
    }

    const auto chrom_has_colon = region.mChromName.find(':') != std::string::npos;
    auto curr_window_start = static_cast<i64>(region.mRegionSpan[0].value());
    const auto max_window_pos = static_cast<i64>(region.mRegionSpan[1].value());

    while ((curr_window_start + window_len) <= max_window_pos) {
      const i64 curr_window_end = curr_window_start + window_len;
      rspec = chrom_has_colon ? fmt::format("{{{}}}:{}-{}", region.mChromName, curr_window_start, curr_window_end)
                              : fmt::format("{}:{}-{}", region.mChromName, curr_window_start, curr_window_end);

      auto wptr = std::make_shared<Window>(mRefPtr->ParseRegion(rspec.c_str()), chrom, mRefPtr->FastaPath());
      uniq_windows.emplace(std::move(wptr));
      curr_window_start += step_size;
    }
  }

  static const auto WindowPtrComparator = [](const WindowPtr &first, const WindowPtr &second) -> bool {
    // NOLINTBEGIN(readability-braces-around-statements)
    if (first->ChromIndex() != second->ChromIndex()) return first->ChromIndex() < second->ChromIndex();
    if (first->StartPos1() != second->StartPos1()) return first->StartPos1() < second->StartPos1();
    return first->EndPos1() < second->EndPos1();
    // NOLINTEND(readability-braces-around-statements)
  };

  usize current_idx = 0;
  std::vector<WindowPtr> results(uniq_windows.cbegin(), uniq_windows.cend());
  std::ranges::sort(results, WindowPtrComparator);
  std::ranges::for_each(results, [&current_idx](WindowPtr &window_ptr) {
    window_ptr->SetGenomeIndex(current_idx);
    current_idx++;
  });

  return results;
}

void WindowBuilder::PadInputRegion(ParseRegionResult &result) const {
  const auto contig_info = mRefPtr->FindChromByName(result.mChromName);
  if (!contig_info.ok()) {
    const auto msg = fmt::format("No chromosome named {} found in reference", result.mChromName);
    throw std::runtime_error(msg);
  }

  const auto contig_max_len = contig_info.value().Length();
  const auto curr_start = result.mRegionSpan[0].value_or(1);
  const auto curr_end = result.mRegionSpan[1].value_or(contig_max_len);

  const auto start_underflows = curr_start <= mParams.mRegionPadding;
  const auto end_overflows = (curr_end > contig_max_len) || ((contig_max_len - curr_end) <= mParams.mRegionPadding);

  start_underflows ? result.mRegionSpan[0] = 1 : result.mRegionSpan[0] = (curr_start - mParams.mRegionPadding);
  end_overflows ? result.mRegionSpan[1] = contig_max_len : result.mRegionSpan[1] = (curr_end + mParams.mRegionPadding);

  if (result.Length() < mParams.mWindowLength) {
    const u64 diff = std::abs(static_cast<i64>(result.Length()) - static_cast<i64>(mParams.mWindowLength) - 1);
    const auto curr_left = result.mRegionSpan[0].value();
    const auto curr_right = result.mRegionSpan[1].value();
    const auto left_new_val = (diff / 2) > curr_left ? curr_left - 1 : curr_left - (diff / 2);
    const auto left_flank = curr_left - left_new_val;
    const auto goes_overmax = curr_right + (diff - left_flank) > contig_max_len;
    result.mRegionSpan[0] = curr_left - left_flank;
    result.mRegionSpan[1] = goes_overmax ? contig_max_len : curr_right + (diff - left_flank);
  }
}

}  // namespace lancet::core

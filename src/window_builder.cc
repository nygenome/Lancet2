#include "lancet2/window_builder.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>

#include "absl/strings/numbers.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "lancet2/log_macros.h"
#include "spdlog/spdlog.h"

namespace lancet2 {
WindowBuilder::WindowBuilder(const std::filesystem::path &ref, u32 region_padding, u32 window_length,
                             u32 pct_window_overlap)
    : refRdr(ref), regionPadding(region_padding), windowLength(window_length), pctWindowOverlap(pct_window_overlap) {}

auto WindowBuilder::AddRegion(const std::string &region_str) -> absl::Status {
  const auto result = ParseRegion(region_str);
  if (!result.ok()) return result.status();
  inputRegions.emplace(result.value().ToRegionString());
  return absl::OkStatus();
}

auto WindowBuilder::AddRegionsBatch(absl::Span<const std::string> region_strs) -> absl::Status {
  for (const auto &regStr : region_strs) {
    auto result = AddRegion(regStr);
    if (!result.ok()) return result;
  }
  return absl::OkStatus();
}

auto WindowBuilder::AddRegionsFromBedFile(const std::filesystem::path &bed) -> absl::Status {
  const auto results = ParseBed(bed);
  if (!results.ok()) return results.status();
  for (const auto &reg : results.value()) {
    inputRegions.insert(reg.ToRegionString());
  }
  return absl::OkStatus();
}

void WindowBuilder::AddAllRefRegions() {
  const auto ctgInfos = refRdr.GetContigs();
  inputRegions.clear();  // Need to clear to ensure we don't end up with duplicate windows at the end
  for (const auto &ctg : ctgInfos) {
    inputRegions.emplace(absl::StrFormat("%s:%d-%d", ctg.contigName, 0, ctg.contigLen));
  }
}

auto WindowBuilder::BuildWindows() const -> absl::StatusOr<std::vector<WindowPtr>> {
  if (IsEmpty()) return absl::FailedPreconditionError("no input regions provided to build windows");

  const auto contigIndexMap = refRdr.GetContigIndexMap();
  absl::flat_hash_set<WindowPtr> uniqPaddedWindowPtrs;
  const auto stepSize = StepSize(pctWindowOverlap, windowLength);
  const auto contigIdxMap = refRdr.GetContigIndexMap();

  for (const auto &regStr : inputRegions) {
    const auto parseResult = ParseRegion(regStr);
    if (!parseResult.ok()) {
      LOG_WARN("%s – Skipping region", parseResult.status().message());
      continue;
    }

    const auto &rawReg = parseResult.value();
    if (!contigIdxMap.contains(rawReg.GetChromName())) {
      throw std::invalid_argument(absl::StrFormat("contig %s is not present in reference", rawReg.GetChromName()));
    }

    const auto paddedResult = PadWindow(rawReg);
    if (!paddedResult.ok()) return paddedResult.status();
    const auto &finalRegion = paddedResult.value();

    if (finalRegion.GetLength() <= windowLength) {
      uniqPaddedWindowPtrs.emplace(std::make_shared<RefWindow>(finalRegion));
      continue;
    }

    i64 currWindowStart = finalRegion.GetStartPos0();
    const auto maxWindowPos = finalRegion.GetEndPos0();

    while (currWindowStart < maxWindowPos) {
      const i64 currWindowEnd = currWindowStart + windowLength;
      auto wPtr = std::make_shared<RefWindow>();
      wPtr->SetChromName(finalRegion.GetChromName());
      wPtr->SetStartPos0(currWindowStart);
      wPtr->SetEndPos0(currWindowEnd);
      uniqPaddedWindowPtrs.emplace(std::move(wPtr));
      currWindowStart += stepSize;
    }
  }

  static const auto Comparator = [&contigIdxMap](const RefWindow &r1, const RefWindow &r2) -> bool {
    if (r1.GetChromName() != r2.GetChromName()) {
      return contigIdxMap.at(r1.GetChromName()) < contigIdxMap.at(r2.GetChromName());
    }
    if (r1.GetStartPos0() != r2.GetStartPos0()) return r1.GetStartPos0() < r2.GetStartPos0();
    return r1.GetEndPos0() < r2.GetEndPos0();
  };

  std::vector<WindowPtr> results(uniqPaddedWindowPtrs.cbegin(), uniqPaddedWindowPtrs.cend());
  std::sort(results.begin(), results.end(),
            [](const WindowPtr &r1, const WindowPtr &r2) -> bool { return Comparator(*r1, *r2); });

  std::for_each(results.begin(), results.end(), [](WindowPtr &w) -> void {
    static usize currWindowIdx = 0;
    w->SetWindowIndex(currWindowIdx);
    currWindowIdx++;
  });

  return results;
}

auto WindowBuilder::StepSize(u32 pct_overlap, u32 window_length) -> i64 {
  const auto rawVal = (static_cast<double>(100 - pct_overlap) / 100.0) * static_cast<double>(window_length);
  // round to ensure that steps always move in multiples of 100
  return static_cast<i64>(std::round(rawVal / 100.0) * 100.0);
}

auto WindowBuilder::ParseRegion(std::string_view region_str) -> absl::StatusOr<RefWindow> {
  std::vector<std::string> tokens = absl::StrSplit(region_str, absl::ByAnyChar(":-"));

  if (tokens.empty() || tokens.size() > 3) {
    const auto errMsg = absl::StrFormat("invalid samtools region string: %s", region_str);
    return absl::InvalidArgumentError(errMsg);
  }

  i64 winStart = 0;
  i64 winEnd = INT64_MAX;

  // NOTE: samtools region strings have 1-based start and end
  if (tokens.size() >= 2) {
    const auto tmp = std::strtoull(tokens[1].c_str(), nullptr, 10);
    winStart = static_cast<i64>(tmp) - 1;
    if (winStart < 0) winStart = 0;
  }

  if (tokens.size() == 3) {
    const auto tmp = std::strtoull(tokens[2].c_str(), nullptr, 10);
    winEnd = static_cast<i64>(tmp) - 1;
  }

  RefWindow w;
  w.SetChromName(tokens[0]);
  w.SetStartPos0(winStart);
  w.SetEndPos0(winEnd);
  return std::move(w);
}

auto WindowBuilder::ParseBed(const std::filesystem::path &bed) -> absl::StatusOr<std::vector<RefWindow>> {
  std::ifstream bedFh(bed, std::ios_base::in);
  std::string line;
  auto lineNum = 0;

  i64 winStart = -1;
  i64 winEnd = -1;
  std::vector<RefWindow> results;

  while (std::getline(bedFh, line)) {
    lineNum++;
    std::vector<std::string> tokens = absl::StrSplit(line, absl::ByChar('\t'), absl::SkipEmpty());

    if (tokens.size() != 3) {
      const auto errMsg = absl::StrFormat("invalid bed line with %d columns at line num %d", tokens.size(), lineNum);
      return absl::InvalidArgumentError(errMsg);
    }

    if (!absl::SimpleAtoi(tokens[1], &winStart) || !absl::SimpleAtoi(tokens[2], &winEnd)) {
      const auto errMsg = absl::StrFormat("could not parse bed line: %s", line);
      return absl::InternalError(errMsg);
    }

    // NOTE: bed file has 0-based start and end
    RefWindow w;
    w.SetChromName(tokens[0]);
    w.SetStartPos0(winStart);
    w.SetEndPos0(winEnd);
    results.emplace_back(std::move(w));
  }

  return std::move(results);
}

auto WindowBuilder::PadWindow(const RefWindow &w) const -> absl::StatusOr<RefWindow> {
  const auto ctgMaxLen = refRdr.GetContigLength(w.GetChromName());
  if (!ctgMaxLen.ok()) return ctgMaxLen.status();

  const auto currMax = static_cast<i64>(ctgMaxLen.value());
  const auto currStart = w.GetStartPos0();
  const auto currEnd = w.GetEndPos0();

  const auto startUnderflows = currStart < regionPadding;
  const auto endOverflows = (currEnd >= currMax) || ((currMax - currEnd) < regionPadding);

  RefWindow result(w);
  startUnderflows ? result.SetStartPos0(0) : result.SetStartPos0(currStart - regionPadding);
  endOverflows ? result.SetEndPos0(currMax) : result.SetEndPos0(currEnd + regionPadding);
  return std::move(result);
}

auto BuildWindows(const absl::flat_hash_map<std::string, i64> &contig_ids, const CliParams &params)
    -> std::vector<WindowPtr> {
  WindowBuilder wb(params.referencePath, params.regionPadLength, params.windowLength, params.pctOverlap);
  for (const auto &region : params.inRegions) {
    const auto result = wb.AddRegion(region);
    if (!result.ok()) {
      LOG_ERROR(result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  if (!params.bedFilePath.empty()) {
    const auto result = wb.AddRegionsFromBedFile(params.bedFilePath);
    if (!result.ok()) {
      LOG_ERROR(result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  if (wb.IsEmpty()) {
    LOG_INFO("No input regions provided to build windows. Using contigs in fasta as input regions");
    wb.AddAllRefRegions();
  }

  LOG_INFO("Building windows using {} input regions", wb.GetSize());
  const auto windows = wb.BuildWindows();
  if (!windows.ok()) {
    LOG_ERROR(windows.status().message());
    std::exit(EXIT_FAILURE);
  }

  return windows.value();
}
}  // namespace lancet2

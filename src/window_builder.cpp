#include "lancet/window_builder.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>

#include "absl/strings/numbers.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "spdlog/spdlog.h"

namespace lancet {
WindowBuilder::WindowBuilder(const std::filesystem::path &ref, std::uint32_t region_padding,
                             std::uint32_t window_length, std::uint32_t pct_window_overlap)
    : refRdr(ref), regionPadding(region_padding), windowLength(window_length), pctWindowOverlap(pct_window_overlap) {}

auto WindowBuilder::AddSamtoolsRegion(const std::string &region_str) -> absl::Status {
  const auto result = ParseRegion(region_str);
  if (!result.ok()) return result.status();
  inputRegions.emplace_back(result.ValueOrDie());
  return absl::OkStatus();
}

auto WindowBuilder::AddBedFileRegions(const std::filesystem::path &bed) -> absl::Status {
  const auto results = ParseBed(bed);
  if (!results.ok()) return results.status();
  const auto &parsedRegions = results.ValueOrDie();
  inputRegions.insert(inputRegions.end(), parsedRegions.begin(), parsedRegions.end());
  return absl::OkStatus();
}

void WindowBuilder::AddAllRefRegions() {
  const auto ctgInfos = refRdr.ContigsInfo();
  for (const auto &ctg : ctgInfos) {
    RefWindow w;
    w.SetChromosome(ctg.contigName);
    w.SetStartPosition0(0);
    w.SetEndPosition0(ctg.contigLen);
    inputRegions.emplace_back(std::move(w));
  }
}

auto WindowBuilder::BuildWindows(const absl::flat_hash_map<std::string, std::int64_t> &contig_ids) const
    -> StatusOr<std::vector<WindowPtr>> {
  if (IsEmpty()) return absl::FailedPreconditionError("no input regions provided to build windows");

  std::vector<WindowPtr> results;
  const auto stepSize = StepSize(pctWindowOverlap, windowLength);

  for (const auto &inRegion : inputRegions) {
    if (!contig_ids.contains(inRegion.Chromosome())) {
      throw std::invalid_argument(absl::StrFormat("contig %s is not present in reference", inRegion.Chromosome()));
    }

    const auto paddedResult = PadWindow(inRegion);
    if (!paddedResult.ok()) return paddedResult.status();
    const auto &region = paddedResult.ValueOrDie();

    if (region.Length() <= windowLength) {
      results.emplace_back(std::make_shared<RefWindow>(region));
      continue;
    }

    std::int64_t currWindowStart = region.StartPosition0();
    const auto maxWindowPos = inRegion.EndPosition0();

    while (currWindowStart < maxWindowPos) {
      const std::int64_t currWindowEnd = currWindowStart + windowLength;
      results.emplace_back(std::make_shared<RefWindow>());
      results.back()->SetChromosome(region.Chromosome());
      results.back()->SetStartPosition0(currWindowStart);
      results.back()->SetEndPosition0(currWindowEnd);
      currWindowStart += stepSize;
    }
  }

  static const auto Comparator = [&contig_ids](const RefWindow &r1, const RefWindow &r2) -> bool {
    if (r1.Chromosome() != r2.Chromosome()) return contig_ids.at(r1.Chromosome()) < contig_ids.at(r2.Chromosome());
    if (r1.StartPosition0() != r2.StartPosition0()) return r1.StartPosition0() < r2.StartPosition0();
    return r1.EndPosition0() < r2.EndPosition0();
  };

  std::sort(results.begin(), results.end(),
            [](const WindowPtr &r1, const WindowPtr &r2) -> bool { return Comparator(*r1, *r2); });

  std::for_each(results.begin(), results.end(), [](WindowPtr &w) -> void {
    static std::size_t currWindowIdx = 0;
    w->SetWindowIndex(currWindowIdx);
    currWindowIdx++;
  });

  return std::move(results);
}

auto WindowBuilder::StepSize(std::uint32_t pct_overlap, std::uint32_t window_length) -> std::int64_t {
  const auto rawVal = (static_cast<double>(100 - pct_overlap) / 100.0) * static_cast<double>(window_length);
  // round to ensure that steps always move in multiples of 100
  return static_cast<std::int64_t>(std::round(rawVal / 100.0) * 100.0);
}

auto WindowBuilder::ParseRegion(std::string_view region_str) -> StatusOr<RefWindow> {
  std::vector<std::string> tokens = absl::StrSplit(region_str, absl::ByAnyChar(":-"));

  if (tokens.empty() || tokens.size() > 3) {
    const auto errMsg = absl::StrFormat("invalid samtools region string: %s", region_str);
    return absl::InvalidArgumentError(errMsg);
  }

  std::int64_t winStart = 0;
  std::int64_t winEnd = INT64_MAX;

  // NOTE: samtools region strings have 1-based start and end
  if (tokens.size() >= 2) {
    const auto tmp = std::strtoull(tokens[1].c_str(), nullptr, 10);
    winStart = static_cast<std::int64_t>(tmp) - 1;
    if (winStart < 0) winStart = 0;
  }

  if (tokens.size() == 3) {
    const auto tmp = std::strtoull(tokens[2].c_str(), nullptr, 10);
    winEnd = static_cast<std::int64_t>(tmp) - 1;
  }

  RefWindow w;
  w.SetChromosome(tokens[0]);
  w.SetStartPosition0(winStart);
  w.SetEndPosition0(winEnd);
  return std::move(w);
}

auto WindowBuilder::ParseBed(const std::filesystem::path &bed) -> StatusOr<std::vector<RefWindow>> {
  std::ifstream bedFh(bed, std::ios_base::in);
  std::string line;
  auto lineNum = 0;

  std::int64_t winStart = -1;
  std::int64_t winEnd = -1;
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
    w.SetChromosome(tokens[0]);
    w.SetStartPosition0(winStart);
    w.SetEndPosition0(winEnd);
    results.emplace_back(std::move(w));
  }

  return std::move(results);
}

auto WindowBuilder::PadWindow(const RefWindow &w) const -> StatusOr<RefWindow> {
  const auto ctgMaxLen = refRdr.ContigLength(w.Chromosome());
  if (!ctgMaxLen.ok()) return ctgMaxLen.status();

  const auto currMax = static_cast<std::int64_t>(ctgMaxLen.ValueOrDie());
  const auto currStart = w.StartPosition0();
  const auto currEnd = w.EndPosition0();

  const auto startUnderflows = currStart < regionPadding;
  const auto endOverflows = (currEnd >= currMax) || ((currMax - currEnd) < regionPadding);

  RefWindow result(w);
  startUnderflows ? result.SetStartPosition0(0) : result.SetStartPosition0(currStart - regionPadding);
  endOverflows ? result.SetEndPosition0(currMax) : result.SetEndPosition0(currEnd + regionPadding);
  return std::move(result);
}

auto BuildWindows(const absl::flat_hash_map<std::string, std::int64_t> &contig_ids, const CliParams &params)
    -> std::vector<WindowPtr> {
  WindowBuilder wb(params.referencePath, params.regionPadLength, params.windowLength, params.pctOverlap);
  const std::vector<std::string> inputRegionStrs = absl::StrSplit(params.inRegionsList, ';');
  for (const auto &region : inputRegionStrs) {
    const auto result = wb.AddSamtoolsRegion(region);
    if (!result.ok()) {
      spdlog::error(result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  if (!params.bedFilePath.empty()) {
    const auto result = wb.AddBedFileRegions(params.bedFilePath);
    if (!result.ok()) {
      spdlog::error(result.message());
      std::exit(EXIT_FAILURE);
    }
  }

  if (wb.IsEmpty()) {
    spdlog::info("No input regions provided to build windows. Using reference contigs in fasta as input");
    wb.AddAllRefRegions();
  }

  spdlog::info("Building windows using {} reference contigs as input regions", wb.Size());
  const auto windows = wb.BuildWindows(contig_ids);
  if (!windows.ok()) {
    spdlog::error(windows.status().message());
    std::exit(EXIT_FAILURE);
  }

  return windows.ValueOrDie();
}
}  // namespace lancet

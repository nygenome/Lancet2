#include "lancet/hts/extractor.h"

#include <algorithm>
#include <filesystem>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

extern "C" {
#include "htslib/cram.h"
#include "htslib/hts.h"
#include "htslib/hts_expr.h"
#include "htslib/hts_log.h"
#include "htslib/sam.h"
}

#include "lancet/base/logging.h"
#include "lancet/base/types.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/reference.h"

#include "absl/container/flat_hash_map.h"
#include "absl/strings/match.h"
#include "absl/strings/str_split.h"
#include "absl/strings/strip.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/core.h"
#include "spdlog/fmt/bundled/format.h"

namespace lancet::hts {

Extractor::Extractor(std::filesystem::path aln_file, Reference const& ref,
                     Alignment::Fields const fields, absl::Span<std::string const> tags,
                     bool const skip_ref_contig_check)
    : mFieldsNeeded(fields),
      mBamCramPath(std::move(aln_file)),
      mTagsNeeded(tags.cbegin(), tags.cend()) {
  hts_set_log_level(HTS_LOG_ERROR);
  auto const bc_path = mBamCramPath.string();

  mFilePtr = InitHtsFile(bc_path.c_str());
  EnsureValidBamOrCram(mFilePtr.get(), bc_path);
  mHdrPtr = InitSamHdr(mFilePtr.get(), bc_path);

  if (!skip_ref_contig_check) HeaderContigsCheck(mHdrPtr.get(), ref);

  SetDefaultHtsOpts(mFilePtr.get(), ref, bc_path);
  mIdxPtr = InitHtsIdx(mFilePtr.get(), bc_path);
  mItrPtr = InitHtsItr(mIdxPtr.get(), bc_path);
  mAlnPtr = InitSamAln(bc_path);
  SetCramRequiredFields(mFieldsNeeded);
}

void Extractor::SetFilterExpression(std::string const& expr) {
  mFiltrPtr.reset(hts_filter_init(expr.c_str()));
  if (mFiltrPtr == nullptr) {
    auto const err_msg = fmt::format("Invalid hts_filter expression: {}", expr);
    throw std::invalid_argument(err_msg);
  }
}

void Extractor::SetRegionToExtract(std::string const& region_spec) {
  mItrPtr.reset(sam_itr_querys(mIdxPtr.get(), mHdrPtr.get(), region_spec.c_str()));
  if (mItrPtr == nullptr) {
    auto const err_msg = fmt::format("Could not set BAM/CRAM iterator for region: {}", region_spec);
    throw std::runtime_error(err_msg);
  }
}

void Extractor::SetRegionToExtract(Reference::Region const& region) {
  SetRegionToExtract(region.ToSamtoolsRegion());
}

void Extractor::SetRegionBatchToExtract(absl::Span<std::string> region_specs) {
  std::vector<char*> regarray;
  regarray.reserve(region_specs.size() + 1);
  std::ranges::transform(region_specs, std::back_inserter(regarray),
                         [](std::string& region) -> char* { return region.data(); });
  regarray.push_back(nullptr);

  mItrPtr.reset(
      sam_itr_regarray(mIdxPtr.get(), mHdrPtr.get(), regarray.data(), region_specs.size()));
  if (mItrPtr == nullptr) {
    throw std::runtime_error("Could not set BAM/CRAM iterator for provided regions");
  }
}

void Extractor::SetRegionBatchToExtract(absl::Span<Reference::Region const> regions) {
  std::vector<std::string> regspecs;
  regspecs.reserve(regions.size());
  std::ranges::transform(regions, std::back_inserter(regspecs),
                         std::mem_fn(&Reference::Region::ToSamtoolsRegion));
  SetRegionBatchToExtract(absl::MakeSpan(regspecs));
}

void Extractor::SetNumThreads(int const nthreads) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
  if (mFilePtr && nthreads > 1 && hts_set_opt(mFilePtr.get(), HTS_OPT_NTHREADS, nthreads) != 0) {
    throw std::runtime_error("Could not set threads");
  }
}

auto Extractor::begin() -> Iterator {
  auto result = Iterator();
  result.mRawFilePtr = mFilePtr.get();
  result.mRawHdrPtr = mHdrPtr.get();
  result.mRawItrPtr = mItrPtr.get();
  result.mRawAlnPtr = mAlnPtr.get();
  result.mRawFiltrPtr = mFiltrPtr.get();
  result.FetchNextAlignment();
  return result;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
auto Extractor::end() -> Iterator {
  return {};
}

auto Extractor::ChromName(i32 chrom_index) const -> std::string {
  auto const* result = sam_hdr_tid2name(mHdrPtr.get(), chrom_index);
  if (result == nullptr) {
    auto const err_msg =
        fmt::format("Reference mIdx {} is not found in BAM/CRAM header", chrom_index);
    throw std::invalid_argument(err_msg);
  }
  return {result};
}

void Extractor::SetCramRequiredFields(Alignment::Fields fields) {
  if (mFilePtr->format.format == cram && fields != Alignment::Fields::AUX_RGAUX) {
    // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-type-vararg)
    cram_set_option(mFilePtr->fp.cram, CRAM_OPT_REQUIRED_FIELDS, fields);
    cram_set_option(mFilePtr->fp.cram, CRAM_OPT_DECODE_MD, 0);
    // NOLINTEND(cppcoreguidelines-pro-type-union-access,cppcoreguidelines-pro-type-vararg)
  }
}

auto Extractor::InitHtsFile(char const* file_path) -> HtsFile {
  auto file_ptr = HtsFile(hts_open(file_path, "r"));
  if (file_ptr == nullptr) {
    auto const err_msg = fmt::format("Could not open alignment file: {}", file_path);
    throw std::runtime_error(err_msg);
  }

  return file_ptr;
}

auto Extractor::InitSamHdr(htsFile* raw_fp, std::string_view aln_path) -> SamHdr {
  auto hdr_ptr = SamHdr(sam_hdr_read(raw_fp));
  if (hdr_ptr == nullptr) {
    auto const err_msg = fmt::format("Cannot read header from BAM/CRAM file: {}", aln_path);
    throw std::runtime_error(err_msg);
  }

  return hdr_ptr;
}

auto Extractor::InitHtsIdx(htsFile* raw_fp, std::string const& aln_path) -> HtsIdx {
  // Try loading alternative mIdx before failing
  auto idx_ptr = HtsIdx(sam_index_load(raw_fp, aln_path.c_str()));
  if (idx_ptr == nullptr) {
    usize const dot_position = aln_path.rfind('.', std::string::npos);
    if (dot_position != 0 && dot_position != std::string::npos) {
      auto const* idx_extension = raw_fp->format.format == cram ? "crai" : "bai";
      auto const alt_idx_path = aln_path.substr(0, dot_position) + idx_extension;
      idx_ptr.reset(sam_index_load2(raw_fp, aln_path.c_str(), alt_idx_path.c_str()));
    }
  }

  if (idx_ptr == nullptr) {
    auto const err_msg = fmt::format("Could not load mIdx for BAM/CRAM: {}", aln_path);
    throw std::runtime_error(err_msg);
  }

  return idx_ptr;
}

auto Extractor::InitHtsItr(hts_idx_t* raw_idx, std::string_view aln_path) -> HtsItr {
  auto itr_ptr = HtsItr(sam_itr_queryi(raw_idx, HTS_IDX_START, 0, 0));
  if (itr_ptr == nullptr) {
    auto const err_msg =
        fmt::format("Could not set BAM/CRAM iterator to start of file: {}", aln_path);
    throw std::runtime_error(err_msg);
  }

  return itr_ptr;
}

auto Extractor::InitSamAln(std::string_view aln_path) -> SamAln {
  auto aln_ptr = SamAln(bam_init1());
  if (aln_ptr == nullptr) {
    auto const err_msg =
        fmt::format("Could not init alignment for BAM/CRAM extractor: {}", aln_path);
    throw std::runtime_error(err_msg);
  }

  return aln_ptr;
}

auto Extractor::ParseSampleName(sam_hdr_t* raw_hdr, std::string_view aln_path) -> std::string {
  static constexpr auto READ_GROUP_PREDICATE = [](std::string_view header_line) -> bool {
    return absl::StartsWith(header_line, "@RG");
  };

  std::vector<std::string_view> const rg_lines =
      absl::StrSplit(raw_hdr->text, '\n', READ_GROUP_PREDICATE);
  absl::flat_hash_map<std::string_view, std::string_view> rg_tags;
  std::string result;

  for (auto const& rgl : rg_lines) {
    if (!absl::StrContains(rgl, "SM")) {
      continue;
    }

    rg_tags.erase(rg_tags.begin(), rg_tags.end());
    rg_tags = absl::StrSplit(absl::StripPrefix(rgl, "@RG\t"), absl::ByAnyChar("\t:"));
    auto const curr_sample = rg_tags.at("SM");

    if (result.empty()) {
      result = curr_sample;
      continue;
    }

    if (curr_sample != result) {
      auto const err_msg =
          fmt::format("Multiple samples in @RG header lines of BAM/CRAM: {}", aln_path);
      throw std::invalid_argument(err_msg);
    }
  }

  return result;
}

void Extractor::EnsureValidBamOrCram(htsFile* raw_fp, std::string_view aln_path) {
  if (raw_fp->format.category != sequence_data) {
    auto const err_msg = fmt::format("Cannot read alignment from non-mDfltSeq data: {}", aln_path);
    throw std::invalid_argument(err_msg);
  }

  if (raw_fp->format.format != bam && raw_fp->format.format != cram) {
    auto* format_description = hts_format_description(&raw_fp->format);
    auto const err_msg =
        fmt::format("Got unexpected alignment file with format: {}", format_description);
    std::free(format_description);  // NOLINT(cppcoreguidelines-no-malloc)
    throw std::invalid_argument(err_msg);
  }

  if (hts_check_EOF(raw_fp) != 1) {
    auto const err_msg =
        fmt::format("BAM/CRAM file {} possibly truncated. Missing EOF block", aln_path);
    throw std::invalid_argument(err_msg);
  }
}

void Extractor::HeaderContigsCheck(sam_hdr_t* raw_hdr, Reference const& ref) {
  auto const chroms = ref.ListChroms();
  if (chroms.size() != static_cast<usize>(sam_hdr_nref(raw_hdr))) {
    LOG_WARN("Number of reference contigs in the BAM/CRAM header don't match the reference FASTA")
  }

  for (auto const& chrom : chroms) {
    auto const hdr_tid = sam_hdr_name2tid(raw_hdr, chrom.Name().c_str());

    if (hdr_tid == -2) {
      throw std::runtime_error("Could not parse BAM/CRAM header");
    }

    if (hdr_tid == -1) {
      auto const err_msg =
          fmt::format("Reference contig {} missing in BAM/CRAM header", chrom.Name());
      throw std::invalid_argument(err_msg);
    }

    auto const hdr_len = static_cast<u64>(sam_hdr_tid2len(raw_hdr, hdr_tid));
    if (hdr_len != chrom.Length()) {
      constexpr auto MISMATCH_FORMAT = "Length mismatch for contig {}. Reference={} BAM/CRAM={}";
      auto const err_msg = fmt::format(MISMATCH_FORMAT, chrom.Name(), chrom.Length(), hdr_len);
      throw std::invalid_argument(err_msg);
    }
  }
}

void Extractor::SetDefaultHtsOpts(htsFile* raw_fp, Reference const& ref,
                                  std::string_view aln_path) {
  static constexpr int DEFAULT_BGZF_CACHE_SIZE = 1 << 24;
  hts_set_cache_size(raw_fp, DEFAULT_BGZF_CACHE_SIZE);

  if (hts_set_fai_filename(raw_fp, ref.FastaPath().c_str()) != 0) {
    auto const err_msg = fmt::format("Could not set reference path to read BAM/CRAM: {}", aln_path);
    throw std::runtime_error(err_msg);
  }
}

}  // namespace lancet::hts

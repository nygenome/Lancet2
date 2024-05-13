#ifndef SRC_LANCET_HTS_EXTRACTOR_H_
#define SRC_LANCET_HTS_EXTRACTOR_H_

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

extern "C" {
#include "htslib/hts.h"
#include "htslib/hts_expr.h"
#include "htslib/sam.h"
}

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"
#include "lancet/base/types.h"
#include "lancet/hts/alignment.h"
#include "lancet/hts/iterator.h"
#include "lancet/hts/reference.h"

namespace lancet::hts {

namespace detail {

struct HtsFileDeleter {
  void operator()(htsFile* file_handle) noexcept { hts_close(file_handle); }
};

struct SamHdrDeleter {
  void operator()(sam_hdr_t* hdr) noexcept { sam_hdr_destroy(hdr); }
};

struct HtsIdxDeleter {
  void operator()(hts_idx_t* idx) noexcept { hts_idx_destroy(idx); }
};

struct HtsItrDeleter {
  void operator()(hts_itr_t* itr) noexcept { hts_itr_destroy(itr); }
};

struct Bam1Deleter {
  void operator()(bam1_t* aln) noexcept { bam_destroy1(aln); }
};

struct HtsFilterDeleter {
  void operator()(hts_filter_t* filter) noexcept { hts_filter_free(filter); }
};

}  // namespace detail

class Extractor {
 public:
  static constexpr auto DEFAULT_FIELDS = Alignment::Fields::AUX_RGAUX;
  Extractor(std::filesystem::path aln_file, const Reference& ref, Alignment::Fields fields = DEFAULT_FIELDS,
            absl::Span<const std::string> tags = {}, bool skip_ref_contig_check = false);

  Extractor() = delete;

  // http://www.htslib.org/doc/samtools.html#FILTER_EXPRESSIONS
  void SetFilterExpression(const std::string& expr);

  void SetRegionToExtract(const std::string& region_spec);
  void SetRegionToExtract(const Reference::Region& region);
  void SetRegionBatchToExtract(absl::Span<std::string> region_specs);
  void SetRegionBatchToExtract(absl::Span<const Reference::Region> regions);

  void SetNumThreads(int nthreads);

  [[nodiscard]] auto begin() -> Iterator;
  [[nodiscard]] auto end() -> Iterator;

  [[nodiscard]] auto ChromName(i32 chrom_index) const -> std::string;
  [[nodiscard]] auto SampleName() const noexcept -> std::string { return mSampleName; }

 private:
  using HtsFile = std::unique_ptr<htsFile, detail::HtsFileDeleter>;
  using SamHdr = std::unique_ptr<sam_hdr_t, detail::SamHdrDeleter>;
  using HtsIdx = std::unique_ptr<hts_idx_t, detail::HtsIdxDeleter>;
  using HtsItr = std::unique_ptr<hts_itr_t, detail::HtsItrDeleter>;
  using HtsFilt = std::unique_ptr<hts_filter_t, detail::HtsFilterDeleter>;
  using SamAln = std::unique_ptr<bam1_t, detail::Bam1Deleter>;

  HtsFile mFilePtr = nullptr;
  SamHdr mHdrPtr = nullptr;
  HtsIdx mIdxPtr = nullptr;
  HtsItr mItrPtr = nullptr;
  HtsFilt mFiltrPtr = nullptr;
  SamAln mAlnPtr = nullptr;
  std::string mSampleName;
  Alignment::Fields mFieldsNeeded;
  std::filesystem::path mBamCramPath;
  absl::flat_hash_set<std::string> mTagsNeeded;

  void SetCramRequiredFields(Alignment::Fields fields);

  [[nodiscard]] static auto InitHtsFile(const char* file_path) -> HtsFile;
  [[nodiscard]] static auto InitSamHdr(htsFile* raw_fp, std::string_view aln_path) -> SamHdr;
  [[nodiscard]] static auto InitHtsIdx(htsFile* raw_fp, const std::string& aln_path) -> HtsIdx;
  [[nodiscard]] static auto InitHtsItr(hts_idx_t* raw_idx, std::string_view aln_path) -> HtsItr;
  [[nodiscard]] static auto InitSamAln(std::string_view aln_path) -> SamAln;

  [[nodiscard]] static auto ParseSampleName(sam_hdr_t* header_line, std::string_view aln_path) -> std::string;

  // Ensure file is a BAM or CRAM and has a valid EOF block
  static void EnsureValidBamOrCram(htsFile* raw_fp, std::string_view aln_path);

  // Check if all contigs present in reference FASTA match the BAM/CRAM mDfltSeq headers
  static void HeaderContigsCheck(sam_hdr_t* raw_hdr, const Reference& ref);

  static void SetDefaultHtsOpts(htsFile* raw_fp, const Reference& ref, std::string_view aln_path);
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_EXTRACTOR_H_

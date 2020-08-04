#include "lancet/hts_reader.h"

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/match.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "lancet/statusor.h"

namespace lancet {
struct BamHdrDeleter {
  void operator()(sam_hdr_t* hdr) noexcept {
    if (hdr != nullptr) bam_hdr_destroy(hdr);
  }
};

struct HtsfileDeleter {
  void operator()(htsFile* sf) noexcept {
    if (sf != nullptr) hts_close(sf);
  }
};

struct HtsIdxDeleter {
  void operator()(hts_idx_t* hidx) noexcept {
    if (hidx != nullptr) hts_idx_destroy(hidx);
  }
};

struct HtsItrDeleter {
  void operator()(hts_itr_t* hitr) noexcept {
    if (hitr != nullptr) hts_itr_destroy(hitr);
  }
};

struct Bam1Deleter {
  void operator()(bam1_t* b) noexcept {
    if (b != nullptr) bam_destroy1(b);
  }
};

static inline auto GetAuxPtr(bam1_t* b, const char* tag) -> StatusOr<const std::uint8_t*> {
  const std::uint8_t* auxData = bam_aux_get(b, tag);
  if (auxData == nullptr && errno == ENOENT) {
    return absl::NotFoundError(absl::StrFormat("could not find tag %s in alignment", tag));
  }
  if (auxData == nullptr) {
    return absl::InternalError(absl::StrFormat("could not get tag %s from alignment", tag));
  }

  return auxData;
}

static const inline auto ShouldSkipAlignment = [](bam1_t* b) -> bool {
  return b == nullptr || (b->core.flag & BAM_FSECONDARY) != 0 || (b->core.flag & BAM_FQCFAIL) != 0 ||  // NOLINT
         (b->core.flag & BAM_FDUP) != 0;                                                               // NOLINT
};

class HtsReader::Impl {
 public:
  Impl(const std::filesystem::path& inpath, const std::filesystem::path& ref) : filePath(inpath) {
    Initialize(inpath, ref);
  }
  ~Impl() = default;
  Impl(Impl&&) noexcept = default;
  auto operator=(Impl&&) noexcept -> Impl& = default;

  Impl() = delete;
  Impl(const Impl&) = delete;
  auto operator=(const Impl&) -> Impl& = delete;

  auto SetRegionSpec(const char* region_spec) -> absl::Status {
    itr.reset(sam_itr_querys(idx.get(), hdr.get(), region_spec));
    if (itr == nullptr) {
      const auto errMsg = absl::StrFormat("could not create iterator to %s for %s", region_spec, filePath);
      return absl::InternalError(errMsg);
    }
    return absl::OkStatus();
  }

  auto SetRegions(absl::Span<const GenomicRegion> regions) -> absl::Status {
    std::vector<std::string> regionStrings;
    for (const auto& region : regions) {
      regionStrings.emplace_back(region.ToRegionString());
    }

    auto regarray = BuildRegionsArray(absl::MakeSpan(regionStrings));
    itr.reset(sam_itr_regarray(idx.get(), hdr.get(), regarray.data(), regions.size()));
    if (itr == nullptr) {
      const auto errMsg = absl::StrFormat("could not create multi-region iterator for %s", filePath);
      return absl::InternalError(errMsg);
    }
    return absl::OkStatus();
  }

  [[nodiscard]] auto NextAlignment(HtsAlignment* result, absl::Span<const std::string> fill_tags) -> IteratorState {
    const auto qryResult = sam_itr_next(fp.get(), itr.get(), aln.get());
    if (qryResult == -1) return IteratorState::DONE;
    if (qryResult < -1) return IteratorState::INVALID;
    if (ShouldSkipAlignment(aln.get())) return NextAlignment(result, fill_tags);

    result->Clear();
    result->SetReadName(bam_get_qname(aln.get()));  // NOLINT
    result->SetContig(sam_hdr_tid2name(hdr.get(), aln->core.tid));
    result->SetMateContig(sam_hdr_tid2name(hdr.get(), aln->core.mtid));
    result->SetStartPosition0(aln->core.pos);
    result->SetMateStartPosition0(aln->core.mpos);
    result->SetEndPosition0(bam_endpos(aln.get()));
    result->SetMappingQuality(aln->core.qual);

    const auto queryLength = static_cast<std::size_t>(aln->core.l_qseq);
    std::string sequence(queryLength, 'N');
    std::string quality(queryLength, static_cast<char>(0));

    const auto* seqBases = bam_get_seq(aln.get());  // NOLINT
    for (std::size_t i = 0; i < queryLength; ++i) {
      sequence[i] = seq_nt16_str[bam_seqi(seqBases, i)];  // NOLINT
    }

    const auto* seqQuals = bam_get_qual(aln.get());  // NOLINT
    for (std::size_t i = 0; i < queryLength; ++i) {
      quality[i] = static_cast<char>(seqQuals[i]);  // NOLINT
    }

    result->SetReadSequence(std::move(sequence));
    result->SetReadQuality(std::move(quality));
    result->SetSamFlags(aln->core.flag);

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wcast-align"
#endif
    const auto* rawCigarData = bam_get_cigar(aln.get());  // NOLINT
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

    const auto cigarLength = static_cast<std::size_t>(aln->core.n_cigar);
    AlignmentCigar cigar;
    cigar.reserve(cigarLength);
    for (std::size_t i = 0; i < cigarLength; ++i) {
      const auto op = bam_cigar_opchr(rawCigarData[i]);                     // NOLINT
      cigar.emplace_back(CigarUnit(op, bam_cigar_oplen(rawCigarData[i])));  // NOLINT
    }
    result->SetCigar(std::move(cigar));

    for (const auto& tag : fill_tags) {
      auto auxResult = GetAuxPtr(aln.get(), tag.c_str());
      if (!auxResult.ok()) continue;
      result->SetTagData(tag, auxResult.ValueOrDie());
    }

    return IteratorState::VALID;
  }

  [[nodiscard]] auto SampleNames() const -> std::vector<std::string> {
    absl::flat_hash_set<std::string> samples;
    absl::flat_hash_map<std::string, std::string> rgTags;

    const auto rgPredicate = [](absl::string_view sv) { return absl::StartsWith(sv, "@RG"); };
    const std::vector<std::string> rgLines = absl::StrSplit(hdr->text, '\n', rgPredicate);

    for (const auto& rgLine : rgLines) {
      if (!absl::StrContains(rgLine, "SM")) continue;
      rgTags.erase(rgTags.begin(), rgTags.end());
      rgTags = absl::StrSplit(absl::StripPrefix(rgLine, "@RG\t"), absl::ByAnyChar("\t:"));
      samples.insert(rgTags.at("SM"));
    }

    return std::vector<std::string>(samples.cbegin(), samples.cend());
  }

  [[nodiscard]] auto ContigsInfo() const -> std::vector<ContigInfo> {
    std::vector<ContigInfo> result;

    const auto numContigs = static_cast<std::size_t>(hdr->n_targets);
    result.reserve(numContigs);

    for (std::size_t i = 0; i < numContigs; ++i) {
      const auto contigName = std::string(hdr->target_name[i]);         // NOLINT
      result.emplace_back(ContigInfo{contigName, hdr->target_len[i]});  // NOLINT
    }

    return result;
  }

  void ResetIterator() {
    itr.reset(sam_itr_queryi(idx.get(), HTS_IDX_START, 0, 0));
    if (itr == nullptr) {
      const auto errMsg = absl::StrFormat("could not set BAM/CRAM iterator to start of file %s", filePath);
      throw std::runtime_error(errMsg);
    }
  }

 private:
  std::filesystem::path filePath;
  std::unique_ptr<htsFile, HtsfileDeleter> fp;
  std::unique_ptr<bam_hdr_t, BamHdrDeleter> hdr;
  std::unique_ptr<hts_idx_t, HtsIdxDeleter> idx;
  std::unique_ptr<hts_itr_t, HtsItrDeleter> itr;
  std::unique_ptr<bam1_t, Bam1Deleter> aln{bam_init1()};

  void Initialize(const std::filesystem::path& inpath, const std::filesystem::path& ref) {
    fp = std::unique_ptr<htsFile, HtsfileDeleter>(hts_open(inpath.c_str(), "r"));
    if (fp == nullptr) {
      const auto errMsg = absl::StrFormat("cannot open hts file %s", inpath);
      throw std::runtime_error(errMsg);
    }

    if (fp->format.category != sequence_data) {
      const auto errMsg = absl::StrFormat("cannot extract reads from non-sequence data %s", inpath);
      throw std::invalid_argument(errMsg);
    }

    if (fp->format.format == sam) {
      const auto errMsg = absl::StrFormat("cannot extract reads from SAM file %s", inpath);
      throw std::invalid_argument(errMsg);
    }

    if (hts_check_EOF(fp.get()) != 1) {
      const auto errMsg = absl::StrFormat("BAM/CRAM file %s possibly truncated. Missing EOF magic block", inpath);
      throw std::invalid_argument(errMsg);
    }

    if (fp->format.format == cram && hts_set_fai_filename(fp.get(), ref.c_str()) != 0) {
      const auto errMsg = absl::StrFormat("could not set reference path %s to read cram %s", ref, inpath);
      throw std::runtime_error(errMsg);
    }

    hdr = std::unique_ptr<bam_hdr_t, BamHdrDeleter>(sam_hdr_read(fp.get()));
    if (hdr == nullptr) {
      const auto errMsg = absl::StrFormat("cannot read header for BAM/CRAM %s", inpath);
      throw std::runtime_error(errMsg);
    }

    // Try loading alternative index before failing
    idx = std::unique_ptr<hts_idx_t, HtsIdxDeleter>(sam_index_load(fp.get(), inpath.c_str()));
    if (idx == nullptr) {
      const auto dotPos = inpath.string().rfind('.', std::string::npos);
      if (dotPos != 0 && dotPos != std::string::npos) {
        const auto* idxExt = fp->format.format == cram ? "crai" : "bai";
        const auto altIdxPath = inpath.string().substr(0, dotPos) + idxExt;
        idx.reset(sam_index_load2(fp.get(), inpath.c_str(), altIdxPath.c_str()));
      }
    }

    if (idx == nullptr) {
      const auto errMsg = absl::StrFormat("could not load index for BAM/CRAM %s", inpath);
      throw std::runtime_error(errMsg);
    }

    ResetIterator();
  }

  static auto BuildRegionsArray(absl::Span<std::string> region_strings) -> std::vector<char*> {
    std::vector<char*> result;
    result.reserve(region_strings.size() + 1);  // +1 for nullptr terminator
    std::transform(region_strings.begin(), region_strings.end(), std::back_inserter(result),
                   [](std::string& s) -> char* { return s.data(); });
    result.push_back(nullptr);
    return result;
  }
};

HtsReader::HtsReader(const std::filesystem::path& inpath, const std::filesystem::path& ref)
    : pimpl(std::make_unique<Impl>(inpath, ref)) {}

HtsReader::~HtsReader() = default;

auto HtsReader::SetRegion(const std::string& contig) -> absl::Status { return pimpl->SetRegionSpec(contig.c_str()); }
auto HtsReader::SetRegion(const GenomicRegion& region) -> absl::Status {
  const auto regionSpec = region.ToRegionString();
  return pimpl->SetRegionSpec(regionSpec.c_str());
}

auto HtsReader::SetRegions(absl::Span<const GenomicRegion> regions) -> absl::Status {
  return pimpl->SetRegions(regions);
}

auto HtsReader::NextAlignment(HtsAlignment* result, absl::Span<const std::string> fill_tags) -> IteratorState {
  return pimpl->NextAlignment(result, fill_tags);
}

auto HtsReader::SampleNames() const -> std::vector<std::string> { return pimpl->SampleNames(); }
auto HtsReader::ContigsInfo() const -> std::vector<ContigInfo> { return pimpl->ContigsInfo(); }
void HtsReader::ResetIterator() { return pimpl->ResetIterator(); }

auto HasTag(const std::filesystem::path& inpath, const std::filesystem::path& ref, const char* tag,
            int max_alignments_to_read) -> bool {
  HtsReader rdr(inpath, ref);
  int currentAlnCnt = 0;

  HtsAlignment aln;
  while (rdr.NextAlignment(&aln, {tag}) == HtsReader::IteratorState::VALID) {
    if (currentAlnCnt == max_alignments_to_read) break;

    currentAlnCnt++;
    if (aln.HasTag(tag)) return true;
  }

  return false;
}
}  // namespace lancet

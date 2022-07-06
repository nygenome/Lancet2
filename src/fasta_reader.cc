#include "lancet2/fasta_reader.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_format.h"
#include "htslib/faidx.h"

namespace lancet2 {
struct FaidxDeleter {
  void operator()(faidx_t* idx) noexcept {
    if (idx != nullptr) fai_destroy(idx);
  }
};

class FastaReader::Impl {
 public:
  explicit Impl(std::filesystem::path ref) : srcPath(std::move(ref)) { LoadIndex(); }

  [[nodiscard]] auto GetRegionSeq(const GenomicRegion& region) const -> absl::StatusOr<std::string> {
    if (!contigLengths.contains(region.GetChromName())) {
      const auto errMsg = absl::StrFormat("contig %s not found in fasta %s", region.GetChromName(), srcPath);
      return absl::InvalidArgumentError(errMsg);
    }

    const auto chromMaxLen = contigLengths.at(region.GetChromName());
    int parsedSeqLen = 0;
    const auto expectedSeqLen = region.GetLength() == 0 ? chromMaxLen : region.GetLength();
    const auto regionString = region.ToSamtoolsRegion();
    char* rawSeqData = fai_fetch(idx.get(), regionString.c_str(), &parsedSeqLen);

    if (rawSeqData == nullptr || static_cast<usize>(parsedSeqLen) != expectedSeqLen) {
      const auto errMsg = absl::StrFormat("could not fetch sequence for %s from %s", regionString, srcPath);
      return rawSeqData == nullptr ? absl::InternalError(errMsg) : absl::FailedPreconditionError(errMsg);
    }

    const auto rawSeqLen = static_cast<usize>(parsedSeqLen);
    std::string resultSeq(rawSeqLen, 'N');
    for (usize baseIdx = 0; baseIdx < rawSeqLen; ++baseIdx) {
      const auto base = absl::ascii_toupper(static_cast<unsigned char>(rawSeqData[baseIdx]));  // NOLINT
      switch (base) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
          resultSeq[baseIdx] = base;
          break;
        default:
          resultSeq[baseIdx] = 'N';
      }
    }

    free(rawSeqData);  // NOLINT
    return std::move(resultSeq);
  }

  [[nodiscard]] auto GetContigs() const -> std::vector<ContigInfo> {
    std::vector<ContigInfo> result;
    result.reserve(numContigs);
    for (const auto& p : contigLengths) {
      result.emplace_back(ContigInfo{p.first, p.second});
    }
    return result;
  }

  [[nodiscard]] auto GetContigIndexMap() const -> absl::flat_hash_map<std::string, i64> { return contigIds; }

  [[nodiscard]] auto GetContigIndex(std::string_view contig) const -> absl::StatusOr<i64> {
    const auto itr = contigIds.find(contig);
    if (itr != contigIds.end()) return itr->second;
    return absl::NotFoundError(absl::StrFormat("contig %s not found", contig));
  }

  [[nodiscard]] auto GetContigLength(std::string_view contig) const -> absl::StatusOr<i64> {
    const auto itr = contigLengths.find(contig);
    if (itr != contigLengths.end()) return itr->second;
    return absl::NotFoundError(absl::StrFormat("contig %s not found", contig));
  }

 private:
  std::filesystem::path srcPath;
  i64 numContigs = 0;
  std::unique_ptr<faidx_t, FaidxDeleter> idx = nullptr;
  absl::flat_hash_map<std::string, i64> contigIds;
  absl::flat_hash_map<std::string, i64> contigLengths;

  void LoadIndex() {
    idx = std::unique_ptr<faidx_t, FaidxDeleter>(fai_load(srcPath.c_str()));
    if (idx == nullptr) {
      const auto errMsg = absl::StrFormat("could not load fasta index from %s", srcPath);
      throw std::runtime_error(errMsg);
    }

    numContigs = faidx_nseq(idx.get());
    for (i64 i = 0; i < numContigs; ++i) {
      const auto ctgName = std::string(faidx_iseq(idx.get(), i));
      const auto ctgLen = faidx_seq_len(idx.get(), ctgName.c_str());

      contigIds[ctgName] = i;
      contigLengths[ctgName] = ctgLen;
    }
  }
};

FastaReader::FastaReader(const std::filesystem::path& ref) { Open(ref); }
FastaReader::~FastaReader() = default;

void FastaReader::Open(const std::filesystem::path& ref) { pimpl = std::make_unique<Impl>(ref); }

auto FastaReader::GetContigSeq(const std::string& contig) const -> absl::StatusOr<std::string> {
  return pimpl->GetRegionSeq(GenomicRegion(contig));
}

auto FastaReader::GetRegionSeq(const GenomicRegion& region) const -> absl::StatusOr<std::string> {
  return pimpl->GetRegionSeq(region);
}

auto FastaReader::GetContigs() const -> std::vector<ContigInfo> { return pimpl->GetContigs(); }

auto FastaReader::GetContigIndexMap() const -> absl::flat_hash_map<std::string, i64> {
  return pimpl->GetContigIndexMap();
}
auto FastaReader::GetContigIndex(std::string_view contig) const -> absl::StatusOr<i64> {
  return pimpl->GetContigIndex(contig);
}

auto FastaReader::GetContigLength(std::string_view contig) const -> absl::StatusOr<i64> {
  return pimpl->GetContigLength(contig);
}
}  // namespace lancet2

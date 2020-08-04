#include "lancet/fasta_reader.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "absl/container/flat_hash_map.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_format.h"
#include "htslib/faidx.h"

namespace lancet {
struct FaidxDeleter {
  void operator()(faidx_t* idx) noexcept {
    if (idx != nullptr) fai_destroy(idx);
  }
};

class FastaReader::Impl {
 public:
  explicit Impl(std::filesystem::path ref) : srcPath(std::move(ref)) { LoadIndex(); }

  [[nodiscard]] auto RegionSequence(const GenomicRegion& region) const -> StatusOr<std::string> {
    if (!contigLengths.contains(region.Chromosome())) {
      const auto errMsg = absl::StrFormat("contig %s not found in fasta %s", region.Chromosome(), srcPath);
      return absl::InvalidArgumentError(errMsg);
    }

    const auto chromMaxLen = contigLengths.at(region.Chromosome());
    int parsedSeqLen = 0;
    const auto expectedSeqLen = region.Length() == 0 ? chromMaxLen : region.Length();
    const auto regionString = region.ToRegionString();
    char* rawSeqData = fai_fetch(idx.get(), regionString.c_str(), &parsedSeqLen);

    if (rawSeqData == nullptr || parsedSeqLen != expectedSeqLen) {
      const auto errMsg = absl::StrFormat("could not fetch sequence for %s from %s", regionString, srcPath);
      return absl::InternalError(errMsg);
    }

    const auto rawSeqLen = static_cast<std::size_t>(parsedSeqLen);
    std::string resultSeq(rawSeqLen, 'N');
    for (std::size_t baseIdx = 0; baseIdx < rawSeqLen; ++baseIdx) {
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

  [[nodiscard]] auto ContigsInfo() const -> std::vector<ContigInfo> {
    std::vector<ContigInfo> result;
    result.reserve(numContigs);
    for (const auto& p : contigLengths) {
      result.emplace_back(ContigInfo{p.first, p.second});
    }
    return result;
  }

  [[nodiscard]] auto ContigId(std::string_view contig) const -> StatusOr<std::int64_t> {
    const auto itr = contigIds.find(contig);
    if (itr != contigIds.end()) return itr->second;
    return absl::NotFoundError(absl::StrFormat("contig %s not found", contig));
  }

  [[nodiscard]] auto ContigLength(std::string_view contig) const -> StatusOr<std::int64_t> {
    const auto itr = contigLengths.find(contig);
    if (itr != contigLengths.end()) return itr->second;
    return absl::NotFoundError(absl::StrFormat("contig %s not found", contig));
  }

 private:
  std::filesystem::path srcPath;
  std::int64_t numContigs = 0;
  std::unique_ptr<faidx_t, FaidxDeleter> idx = nullptr;
  absl::flat_hash_map<std::string, std::int64_t> contigIds;
  absl::flat_hash_map<std::string, std::int64_t> contigLengths;

  void LoadIndex() {
    idx = std::unique_ptr<faidx_t, FaidxDeleter>(fai_load(srcPath.c_str()));
    if (idx == nullptr) {
      const auto errMsg = absl::StrFormat("could not load fasta index from %s", srcPath);
      throw std::runtime_error(errMsg);
    }

    numContigs = faidx_nseq(idx.get());
    for (std::int64_t i = 0; i < numContigs; ++i) {
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

auto FastaReader::ContigSequence(const std::string& contig) const -> StatusOr<std::string> {
  return pimpl->RegionSequence(GenomicRegion(contig));
}

auto FastaReader::RegionSequence(const GenomicRegion& region) const -> StatusOr<std::string> {
  return pimpl->RegionSequence(region);
}

auto FastaReader::ContigsInfo() const -> std::vector<ContigInfo> { return pimpl->ContigsInfo(); }

auto FastaReader::ContigId(std::string_view contig) const -> StatusOr<std::int64_t> { return pimpl->ContigId(contig); }

auto FastaReader::ContigLength(std::string_view contig) const -> StatusOr<std::int64_t> {
  return pimpl->ContigLength(contig);
}
}  // namespace lancet

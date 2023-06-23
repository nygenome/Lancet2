#ifndef SRC_LANCET_HTS_REFERENCE_H_
#define SRC_LANCET_HTS_REFERENCE_H_

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

extern "C" {
#include "htslib/faidx.h"
}

#include "absl/status/statusor.h"
#include "lancet/base/assert.h"
#include "lancet/base/types.h"

namespace lancet::hts {

namespace detail {

struct FaidxDeleter {
  void operator()(faidx_t* idx) noexcept { fai_destroy(idx); }
};

}  // namespace detail

class Reference {
 public:
  Reference(std::filesystem::path reference);
  Reference() = delete;

  [[nodiscard]] auto FastaPath() const noexcept -> std::filesystem::path { return mFastaPath; }

  class Chrom;

  [[nodiscard]] auto ListChroms() const noexcept -> std::vector<Chrom>;
  [[nodiscard]] auto FindChromByName(std::string_view chrom_name) const noexcept -> absl::StatusOr<Chrom>;
  [[nodiscard]] auto FindChromByIndex(i64 chrom_index) const noexcept -> absl::StatusOr<Chrom>;

  class Region;

  // 1-based inclusive closed position coordinates. Use std::nullopt to skip providing a co-ordinate.
  using OneBasedClosedOptional = std::array<std::optional<u64>, 2>;
  static constexpr auto NULL_INTERVAL = OneBasedClosedOptional{std::nullopt, std::nullopt};

  struct ParseRegionResult {
    // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
    std::string mChromName;
    OneBasedClosedOptional mRegionSpan = NULL_INTERVAL;
    // NOLINTEND(misc-non-private-member-variables-in-classes)

    [[nodiscard]] auto Length() const noexcept -> usize;

    template <typename HashState>
    friend auto AbslHashValue(HashState state, const ParseRegionResult& result) noexcept -> HashState {
      return HashState::combine(std::move(state), result.mChromName, result.mRegionSpan);
    }

    auto operator==(const ParseRegionResult& rhs) const -> bool {
      return mChromName == rhs.mChromName && mRegionSpan == rhs.mRegionSpan;
    }

    auto operator!=(const ParseRegionResult& rhs) const -> bool { return !(rhs == *this); }
  };

  [[nodiscard]] auto ParseRegion(const char* region_spec) const -> ParseRegionResult;

  // Makes a `Reference::Region` after removing any shorthand for start and end co-ordinates.
  // Throws if 0-based co-ordinates given (or) any of the co-ordinates > chromosome mVarLength (or) if end is < start.
  [[nodiscard]] auto MakeRegion(const std::string& chrom_name, const OneBasedClosedOptional& interval) const -> Region;
  [[nodiscard]] auto MakeRegion(const ParseRegionResult& parse_result) const -> Region;
  [[nodiscard]] auto MakeRegion(const char* region_spec) const -> Region;

  template <typename HashState>
  friend auto AbslHashValue(HashState state, const Reference& ref) noexcept -> HashState {
    return HashState::combine(std::move(state), ref.mChroms);
  }

 private:
  using FastaIndex = std::unique_ptr<faidx_t, detail::FaidxDeleter>;

  std::filesystem::path mFastaPath;
  FastaIndex mFastaIndex = nullptr;
  std::vector<Chrom> mChroms;

  using OneBasedClosedInterval = std::array<u64, 2>;
  [[nodiscard]] auto FetchSeq(const std::string& chrom, const OneBasedClosedInterval& full_intvl) const -> std::string;
};

class Reference::Chrom {
 public:
  Chrom() = default;

  [[nodiscard]] auto Name() const noexcept -> std::string { return mName; }
  [[nodiscard]] auto Index() const noexcept -> usize { return mIdx; }
  [[nodiscard]] auto Length() const noexcept -> u64 { return mLength; }

  template <typename HashState>
  friend auto AbslHashValue(HashState state, const Reference::Chrom& chrom) noexcept -> HashState {
    return HashState::combine(std::move(state), chrom.mIdx, chrom.mLength, chrom.mName);
  }

  auto operator==(const Chrom& rhs) const -> bool { return mIdx == rhs.mIdx && mLength == rhs.mLength; }
  auto operator!=(const Chrom& rhs) const -> bool { return mIdx != rhs.mIdx || mLength != rhs.mLength; }
  auto operator<(const Chrom& rhs) const -> bool { return mIdx < rhs.mIdx; }
  auto operator>(const Chrom& rhs) const -> bool { return mIdx > rhs.mIdx; }
  auto operator<=(const Chrom& rhs) const -> bool { return mIdx <= rhs.mIdx; }
  auto operator>=(const Chrom& rhs) const -> bool { return mIdx >= rhs.mIdx; }

 private:
  usize mIdx = -1;
  u64 mLength = 0;
  std::string mName;

  friend class Reference;

  Chrom(i32 chrom_index, std::string_view chrom_name, hts_pos_t chrom_len)
      : mIdx(static_cast<usize>(chrom_index)), mLength(static_cast<u64>(chrom_len)), mName(chrom_name) {}
};

class Reference::Region {
 public:
  // Region specified as: RNAME[:STARTPOS[-ENDPOS]]. Both coordinates are 1-based. i.e. one based closed coordinates.
  // e.g: `chr1`, `chr1:1000-1200`, `chr1:100-` -> short for `chr1:100-${END}`, `chr1:-100` -> short for `chr1:1-100`
  [[nodiscard]] auto ToSamtoolsRegion() const -> std::string;

  [[nodiscard]] auto ChromName() const -> std::string { return mName; }
  [[nodiscard]] auto ChromIndex() const -> usize { return mChromIdx; }
  [[nodiscard]] auto StartPos1() const -> u64 { return mStart1; }
  [[nodiscard]] auto EndPos1() const -> u64 { return mEnd1; }
  [[nodiscard]] auto SeqView() const -> std::string_view { return mSeq; }
  [[nodiscard]] auto SeqData() const -> const char* { return mSeq.c_str(); }

  [[nodiscard]] auto Length() const -> u64 {
    LANCET_ASSERT(mSeq.length() == (mEnd1 - mStart1 + 1))
    return mEnd1 - mStart1 + 1;
  }

  template <typename HashState>
  friend auto AbslHashValue(HashState state, const Reference::Region& reg) noexcept -> HashState {
    return HashState::combine(std::move(state), reg.mChromIdx, reg.mStart1, reg.mEnd1, reg.mName, reg.mSeq);
  }

 private:
  usize mChromIdx = -1;
  u64 mStart1 = 0;
  u64 mEnd1 = 0;
  std::string mName;
  std::string mSeq;

  friend class Reference;

  Region(usize chrom_index, const std::pair<u64, u64>& interval, const char* name, std::string&& seq)
      : mChromIdx(chrom_index), mStart1(interval.first), mEnd1(interval.second), mName(name), mSeq(std::move(seq)) {}
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_REFERENCE_H_

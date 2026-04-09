#ifndef SRC_LANCET_HTS_AUX_TAG_H_
#define SRC_LANCET_HTS_AUX_TAG_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include <cmath>

extern "C" {
#include "htslib/sam.h"
}

#include "lancet/base/types.h"

#include "absl/container/fixed_array.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "spdlog/fmt/bundled/format.h"
#include "spdlog/fmt/bundled/ranges.h"

namespace lancet::hts {

class AuxTag {
 public:
  [[nodiscard]] auto Name() const noexcept -> std::string_view { return {mTagName.data(), 2}; }

  template <typename ResultType>
  [[nodiscard]] auto Value() const noexcept -> absl::StatusOr<ResultType> {
    StaticAssertResultType<ResultType>();
    return GetResultIfAvailable<ResultType>();
  }

  [[nodiscard]] auto ToString() const noexcept -> std::string { return FormatToSamAux(); }

  template <typename HashState>
  friend auto AbslHashValue(HashState state, AuxTag const& aux) noexcept -> HashState {
    return HashState::combine(std::move(state), aux.mIsSigned, aux.mTagName, aux.mCharData,
                              aux.mIntData, aux.mFloatData, aux.mStrData, aux.mArrIntData,
                              aux.mArrFloatData);
  }

  auto operator==(AuxTag const& rhs) const -> bool {
    return mTagName == rhs.mTagName &&
           mCharData == rhs.mCharData &&
           mIntData == rhs.mIntData &&
           mFloatData == rhs.mFloatData &&
           mStrData == rhs.mStrData &&
           mArrIntData == rhs.mArrIntData &&
           mArrFloatData == rhs.mArrFloatData;
  }

  auto operator!=(AuxTag const& rhs) const -> bool { return !(rhs == *this); }

 private:
  static constexpr int MISSING_CHAR = std::numeric_limits<int>::max();
  static constexpr i64 MISSING_INT = std::numeric_limits<i64>::max();
  static constexpr f64 MISSING_FLOAT = std::numeric_limits<f64>::max();
  static constexpr f64 EPSILON = std::numeric_limits<f64>::epsilon();
  using IntVector = absl::FixedArray<i64>;
  using FloatVector = absl::FixedArray<f64>;

  bool mIsSigned = true;
  std::array<char, 2> mTagName = {'\0', '\0'};
  int mCharData = MISSING_CHAR;
  i64 mIntData = MISSING_INT;
  f64 mFloatData = MISSING_FLOAT;
  std::shared_ptr<std::string> mStrData = nullptr;
  std::shared_ptr<IntVector> mArrIntData = nullptr;
  std::shared_ptr<FloatVector> mArrFloatData = nullptr;

  friend class Alignment;

  template <typename T>
  void PopulateArrayData(u8 const* data, usize const arr_len) {
    if constexpr (std::is_same<i64, T>::value) {
      mArrIntData = std::make_shared<IntVector>(arr_len, MISSING_INT);
      for (usize idx = 0; idx < arr_len; ++idx) {
        mArrIntData->at(idx) = bam_auxB2i(data, idx);
      }
    }

    if constexpr (std::is_same<f64, T>::value) {
      mArrFloatData = std::make_shared<FloatVector>(arr_len, MISSING_FLOAT);
      for (usize idx = 0; idx < arr_len; ++idx) {
        mArrFloatData->at(idx) = bam_auxB2f(data, idx);
      }
    }
  }

  explicit AuxTag(u8 const* data) {
    mTagName.at(0) = static_cast<char>(data[-2]);
    mTagName.at(1) = static_cast<char>(data[-1]);

    auto const tag_type = static_cast<char>(*data);
    mIsSigned = (tag_type != 'C' && tag_type != 'S' && tag_type != 'I');

    // NOLINTNEXTLINE(bugprone-switch-missing-default-case)
    switch (tag_type) {
      case 'A':
        // NOLINTNEXTLINE(bugprone-signed-char-misuse,cert-str34-c)
        mCharData = static_cast<int>(bam_aux2A(data));
        break;

      case 'c':
      case 'C':
      case 's':
      case 'S':
      case 'i':
      case 'I':
        mIntData = bam_aux2i(data);
        break;

      case 'f':
      case 'd':
        mFloatData = bam_aux2f(data);
        break;

      case 'Z':
      case 'H':
        mStrData = std::make_shared<std::string>(bam_aux2Z(data));
        break;

      case 'B':
        auto const arr_len = static_cast<usize>(bam_auxB_len(data));
        auto const tag_subtype = static_cast<char>(data[1]);
        mIsSigned = (tag_subtype != 'C' && tag_subtype != 'S' && tag_subtype != 'I');

        // NOLINTNEXTLINE(bugprone-switch-missing-default-case)
        switch (tag_subtype) {
          case 'c':
          case 'C':
          case 's':
          case 'S':
          case 'i':
          case 'I':
            PopulateArrayData<i64>(data, arr_len);
            break;

          case 'f':
          case 'd':
            PopulateArrayData<f64>(data, arr_len);
            break;
        }
        break;
    }
  }

  template <typename ResultType>
  static constexpr void StaticAssertResultType() {
    constexpr auto is_char = std::is_same_v<char, ResultType>;
    constexpr auto is_int_64 = std::is_same_v<i64, ResultType>;
    constexpr auto is_float_64 = std::is_same_v<f64, ResultType>;
    constexpr auto is_string_view = std::is_same_v<std::string_view, ResultType>;
    constexpr auto is_span_int_64 = std::is_same_v<absl::Span<i64 const>, ResultType>;
    constexpr auto is_span_float_64 = std::is_same_v<absl::Span<f64 const>, ResultType>;

    static_assert(
        is_char || is_int_64 || is_float_64 || is_string_view || is_span_int_64 || is_span_float_64,
        "AuxData's Value ResultType must be one of the following types - char, int64_t, double, "
        "std::string_view, absl::Span<const int64_t>, absl::Span<const double>");
  }

  template <typename ResultType>
  [[nodiscard]] static auto MakeInvalidTypeStatus(std::array<char, 2> const tag_name)
      -> absl::Status {
    std::string_view const name_view(tag_name.data(), 2);
    return {absl::StatusCode::kInvalidArgument,
            fmt::format("Tag {} does not have data with type {}", name_view,
                        typeid(ResultType).name())};
  }

  template <typename ResultType>
  [[nodiscard]] auto GetResultIfAvailable() const noexcept -> absl::StatusOr<ResultType> {
    if constexpr (std::is_same<char, ResultType>::value) {
      if (mCharData == MISSING_CHAR)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return static_cast<char>(mCharData);
    }

    if constexpr (std::is_same<i64, ResultType>::value) {
      if (mIntData == MISSING_INT)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return mIntData;
    }

    if constexpr (std::is_same<f64, ResultType>::value) {
      if (std::fabs(mFloatData - MISSING_FLOAT) <= EPSILON)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return mFloatData;
    }

    if constexpr (std::is_same<std::string_view, ResultType>::value) {
      if (mStrData == nullptr)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return std::string_view(*mStrData);
    }

    if constexpr (std::is_same<absl::Span<i64 const>, ResultType>::value) {
      if (mArrIntData == nullptr)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return absl::MakeConstSpan(*mArrIntData);
    }
    if constexpr (std::is_same<absl::Span<f64 const>, ResultType>::value) {
      if (mArrFloatData == nullptr)
        return MakeInvalidTypeStatus<ResultType>(mTagName);
      return absl::MakeConstSpan(*mArrFloatData);
    }
  }

  [[nodiscard]] auto FormatToSamAux() const noexcept -> std::string {
    if (mCharData != MISSING_CHAR) {
      return fmt::format("{}:A:{}", std::string_view(mTagName.data(), 2),
                         static_cast<char>(mCharData));
    }

    if (mIntData != MISSING_INT) {
      char const itype = mIsSigned ? 'i' : 'I';
      return fmt::format("{}:{}:{}", std::string_view(mTagName.data(), 2), itype, mIntData);
    }

    if (mFloatData != MISSING_FLOAT) {
      return fmt::format("{}:f:{}", std::string_view(mTagName.data(), 2), mFloatData);
    }

    if (mStrData != nullptr) {
      return fmt::format("{}:Z:{}", std::string_view(mTagName.data(), 2), *mStrData);
    }

    if (mArrIntData != nullptr) {
      char const subtype = mIsSigned ? 'i' : 'I';
      return fmt::format("{}:B:{}{}", std::string_view(mTagName.data(), 2), subtype,
                         fmt::join(*mArrIntData, ","));
    }

    if (mArrFloatData != nullptr) {
      return fmt::format("{}:B:f{}", std::string_view(mTagName.data(), 2),
                         fmt::join(*mArrFloatData, ","));
    }

    return {};
  }
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_AUX_TAG_H_

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

namespace lancet {
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-const-variable"
#endif
static constexpr char ALIGNMENT_GAP = '-';
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

struct AlignedSequences {
  std::string ref;
  std::string qry;
};

struct AlignedSequencesView {
  std::string_view ref{};
  std::string_view qry{};
};

[[nodiscard]] auto Align(std::string_view ref, std::string_view qry) -> AlignedSequences;

namespace detail {
template <typename T>
class NaiveTwoDimArray {
 public:
  explicit NaiveTwoDimArray(std::size_t nrows, std::size_t ncols)
      : numRows(nrows), numCols(ncols), data(std::make_unique<T[]>(nrows * ncols)) {}  // NOLINT

  NaiveTwoDimArray() = delete;

  auto at(std::size_t row, std::size_t col) -> T& { return data[row * numCols + col]; }
  auto at(std::size_t row, std::size_t col) const -> const T& { return data[row * numCols + col]; }

 private:
  std::size_t numRows;
  std::size_t numCols;
  std::unique_ptr<T[]> data;  // NOLINT
};
}  // namespace detail
}  // namespace lancet

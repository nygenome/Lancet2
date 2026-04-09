#ifndef SRC_LANCET_HTS_ITERATOR_H_
#define SRC_LANCET_HTS_ITERATOR_H_

#include <iterator>

extern "C" {
#include "htslib/hts.h"
#include "htslib/hts_expr.h"
#include "htslib/sam.h"
}

#include "lancet/base/types.h"
#include "lancet/hts/alignment.h"

namespace lancet::hts {

class Iterator {
 public:
  // NOLINTBEGIN(readability-identifier-naming)  // std::iterator_traits mandates snake_case
  using iterator_category = std::input_iterator_tag;
  using value_type = Alignment const;
  using difference_type = std::ptrdiff_t;
  using pointer = Alignment const*;
  using reference = Alignment const&;

  auto operator*() -> reference { return mParsedAln; }
  auto operator==(Iterator const& rhs) const -> bool;
  auto operator!=(Iterator const& rhs) const -> bool;

  auto operator++() -> Iterator&;
  auto operator++(int) -> Iterator;
  // NOLINTEND(readability-identifier-naming)

 private:
  Alignment mParsedAln;
  htsFile* mRawFilePtr = nullptr;
  sam_hdr_t* mRawHdrPtr = nullptr;
  hts_itr_t* mRawItrPtr = nullptr;
  bam1_t* mRawAlnPtr = nullptr;
  hts_filter_t* mRawFiltrPtr = nullptr;

  friend class Extractor;

  Iterator() = default;

  void FetchNextAlignment();
  [[nodiscard]] auto PassesFilter() const -> bool;
};

}  // namespace lancet::hts

#endif  // SRC_LANCET_HTS_ITERATOR_H_

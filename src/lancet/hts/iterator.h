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
  using iterator_category = std::input_iterator_tag;
  using value_type = const Alignment;
  using difference_type = usize;
  using pointer = const Alignment*;
  using reference = const Alignment&;

  auto operator*() -> reference { return mParsedAln; }
  auto operator==(const Iterator& rhs) const -> bool;
  auto operator!=(const Iterator& rhs) const -> bool;

  auto operator++() -> Iterator&;
  auto operator++(int) -> Iterator;

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

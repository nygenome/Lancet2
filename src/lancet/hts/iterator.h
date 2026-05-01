#ifndef SRC_LANCET_HTS_ITERATOR_H_
#define SRC_LANCET_HTS_ITERATOR_H_

#include "lancet/hts/alignment.h"

extern "C" {
#include "htslib/hts.h"
#include "htslib/hts_expr.h"
#include "htslib/sam.h"
}

#include <iterator>
#include <string>
#include <string_view>

#include <cstddef>

namespace lancet::hts {

/// BAM/CRAM record iterator over a region set by Extractor.
///
/// WARNING — Alignment lifetime contract:
/// Alignment objects returned by dereferencing this iterator (`*itr`) are
/// invalidated when the iterator advances (`++itr`). The underlying bam1_t*
/// is reused across steps to avoid per-record allocation. Consequently:
///   - Do NOT store Alignment objects across iterator increments.
///   - Do NOT hold string_views from QnameView() across increments.
///   - Use BuildSequence()/BuildQualities() to deep-copy data
///     that must outlive the current step.
class Iterator {
 public:
  // std::iterator_traits mandates snake_case names for the trait aliases below.
  // NOLINTBEGIN(readability-identifier-naming)
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
  // ── 8B Align ────────────────────────────────────────────────────────────
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

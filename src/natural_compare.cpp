#include "lancet/natural_compare.h"

#include <cctype>
#include <cstring>

namespace lancet {
namespace detail {
static inline auto IsDigit(const char c) -> bool {
  return static_cast<bool>(std::isdigit(static_cast<unsigned char>(c)));
}

static inline auto IsSpace(const char c) -> bool {
  return static_cast<bool>(std::isspace(static_cast<unsigned char>(c)));
}

static auto CompareRight(const char *a, const char *b) -> int {
  int bias = 0;

  /* The longest run of digits wins.  That aside, the greatest
 value wins, but we can't know that it will until we've scanned
 both numbers to know that they have the same magnitude, so we
 remember it in BIAS. */
  for (;; a++, b++) {  // NOLINT
    if (!IsDigit(*a) && !IsDigit(*b)) {
      return bias;
    }
    if (!IsDigit(*a)) {
      return -1;
    }
    if (!IsDigit(*b)) {
      return +1;
    }
    if (*a < *b) {
      if (!static_cast<bool>(bias)) {
        bias = -1;
      }
    } else if (*a > *b) {
      if (!static_cast<bool>(bias)) {
        bias = +1;
      }
    } else if (!static_cast<bool>(*a) && !static_cast<bool>(*b)) {
      return bias;
    }
  }
}

static auto CompareLeft(const char *a, const char *b) -> int {
  /* Compare two left-aligned numbers: the first to have a
     different value wins. */
  for (;; a++, b++) {  // NOLINT
    if (!IsDigit(*a) && !IsDigit(*b)) {
      return 0;
    }
    if (!IsDigit(*a)) {
      return -1;
    }
    if (!IsDigit(*b)) {
      return +1;
    }
    if (*a < *b) {
      return -1;
    }
    if (*a > *b) {
      return +1;
    }
  }
}

static auto StrNatCmp(const char *a, const char *b) -> int {
  int ai = 0;
  int bi = 0;
  char ca = '\0';
  char cb = '\0';
  int fractional = 0;
  int result = 0;

  while (true) {
    ca = a[ai];  // NOLINT
    cb = b[bi];  // NOLINT

    /* skip over leading spaces or zeros */
    while (IsSpace(ca)) {
      ca = a[++ai];  // NOLINT
    }
    while (IsSpace(cb)) {
      cb = b[++bi];  // NOLINT
    }

    /* process run of digits */
    if (IsDigit(ca) && IsDigit(cb)) {
      fractional = static_cast<int>(ca == '0' || cb == '0');

      if (static_cast<bool>(fractional)) {
        if ((result = CompareLeft(a + ai, b + bi)) != 0) {  // NOLINT
          return result;
        }
      } else {
        if ((result = CompareRight(a + ai, b + bi)) != 0) {  // NOLINT
          return result;
        }
      }
    }

    if (!static_cast<bool>(ca) && !static_cast<bool>(cb)) {
      // The strings compare the same.
      // call strcmp to break the tie.
      return std::strcmp(a, b);
    }

    if (ca < cb) {
      return -1;
    }
    if (ca > cb) {
      return +1;
    }

    ++ai;
    ++bi;
  }
}
}  // namespace detail

auto NaturalCompare(const std::string &lhs, const std::string &rhs) -> int {
  return detail::StrNatCmp(lhs.c_str(), rhs.c_str());
}
}  // namespace lancet

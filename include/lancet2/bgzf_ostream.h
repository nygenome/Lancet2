#pragma once

#include <ostream>
#include <string>

#include "lancet2/bgzf_streambuf.h"

namespace lancet2 {
enum class BgzfFormat { UNSPECIFIED, GFF, BED, VCF };

class BgzfOstream : public std::ostream {
 public:
  BgzfOstream() : std::ostream(nullptr) {}
  ~BgzfOstream() override { close(); }

  BgzfOstream(const BgzfOstream &) = delete;
  BgzfOstream(BgzfOstream &&) = delete;
  auto operator=(const BgzfOstream &) -> BgzfOstream & = delete;
  auto operator=(BgzfOstream &&) -> BgzfOstream & = delete;

  auto open(const std::string &path, BgzfFormat ofmt) -> bool;
  auto open(const std::string &path) -> bool { return open(path, BgzfFormat::UNSPECIFIED); }

  void close();

 protected:
  detail::BgzfStreambuf buf;                    // NOLINT
  BgzfFormat outFmt = BgzfFormat::UNSPECIFIED;  // NOLINT

  void build_index();
};
}  // namespace lancet2

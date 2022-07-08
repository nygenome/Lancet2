#include "lancet2/node_label.h"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("Can build and merge adjacent node labels", "[lancet2::NodeLabel]") {
  lancet2::NodeLabel l1(5);
  lancet2::NodeLabel l2(5);

  l1.Push(lancet2::KmerLabel::NORMAL);
  l2.Push(lancet2::KmerLabel::TUMOR);

  CHECK_FALSE(l1.HasLabel(lancet2::KmerLabel::TUMOR));
  l1.MergeBuddy(l2, lancet2::BuddyPosition::FRONT, false, 3);
  CHECK(l1.HasLabel(lancet2::KmerLabel::TUMOR));
}

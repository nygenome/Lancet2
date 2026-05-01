#include "lancet/cbdg/sample_mask.h"

#include "catch_amalgamated.hpp"

using lancet::cbdg::SampleMask;

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SampleMask::SetBit and TestBit", "[lancet][cbdg][SampleMask]") {
  SampleMask mask;

  SECTION("Default mask has no bits set") {
    CHECK_FALSE(mask.TestBit(0));
    CHECK_FALSE(mask.TestBit(1));
    CHECK_FALSE(mask.TestBit(63));
    CHECK_FALSE(mask.TestBit(64));
    CHECK(mask.PopCount() == 0);
  }

  SECTION("Set and test reference bit (bit 0)") {
    mask.SetBit(0);
    CHECK(mask.TestBit(0));
    CHECK_FALSE(mask.TestBit(1));
    CHECK(mask.PopCount() == 1);
  }

  SECTION("Set and test single sample bit") {
    mask.SetBit(5);
    CHECK(mask.TestBit(5));
    CHECK_FALSE(mask.TestBit(4));
    CHECK_FALSE(mask.TestBit(6));
    CHECK(mask.PopCount() == 1);
  }

  SECTION("Set multiple bits within single word") {
    mask.SetBit(0);
    mask.SetBit(1);
    mask.SetBit(10);
    CHECK(mask.TestBit(0));
    CHECK(mask.TestBit(1));
    CHECK(mask.TestBit(10));
    CHECK_FALSE(mask.TestBit(2));
    CHECK(mask.PopCount() == 3);
  }

  SECTION("Setting same bit twice is idempotent") {
    mask.SetBit(7);
    mask.SetBit(7);
    CHECK(mask.TestBit(7));
    CHECK(mask.PopCount() == 1);
  }
}

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SampleMask word boundary behavior", "[lancet][cbdg][SampleMask]") {
  SampleMask mask;

  SECTION("Bit 63 stays in first word") {
    mask.SetBit(63);
    CHECK(mask.TestBit(63));
    CHECK_FALSE(mask.TestBit(62));
    CHECK_FALSE(mask.TestBit(64));
    CHECK(mask.PopCount() == 1);
  }

  SECTION("Bit 64 triggers second word allocation") {
    mask.SetBit(64);
    CHECK(mask.TestBit(64));
    CHECK_FALSE(mask.TestBit(63));
    CHECK_FALSE(mask.TestBit(65));
    CHECK(mask.PopCount() == 1);
  }

  SECTION("Bits spanning two words") {
    mask.SetBit(0);
    mask.SetBit(63);
    mask.SetBit(64);
    mask.SetBit(127);
    CHECK(mask.TestBit(0));
    CHECK(mask.TestBit(63));
    CHECK(mask.TestBit(64));
    CHECK(mask.TestBit(127));
    CHECK_FALSE(mask.TestBit(1));
    CHECK_FALSE(mask.TestBit(128));
    CHECK(mask.PopCount() == 4);
  }
}

// Catch2 SECTION fan-out inflates clang-tidy's cognitive-complexity metric beyond the project
// ceiling.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SampleMask::Merge", "[lancet][cbdg][SampleMask]") {
  SECTION("Merge two non-overlapping masks") {
    SampleMask mask_a;
    SampleMask mask_b;
    mask_a.SetBit(1);
    mask_a.SetBit(3);
    mask_b.SetBit(2);
    mask_b.SetBit(4);

    mask_a.Merge(mask_b);
    CHECK(mask_a.TestBit(1));
    CHECK(mask_a.TestBit(2));
    CHECK(mask_a.TestBit(3));
    CHECK(mask_a.TestBit(4));
    CHECK(mask_a.PopCount() == 4);
  }

  SECTION("Merge overlapping masks is idempotent for shared bits") {
    SampleMask mask_a;
    SampleMask mask_b;
    mask_a.SetBit(1);
    mask_a.SetBit(5);
    mask_b.SetBit(5);
    mask_b.SetBit(10);

    mask_a.Merge(mask_b);
    CHECK(mask_a.TestBit(1));
    CHECK(mask_a.TestBit(5));
    CHECK(mask_a.TestBit(10));
    CHECK(mask_a.PopCount() == 3);
  }

  SECTION("Merge grows target when source has more words") {
    SampleMask mask_a;
    SampleMask mask_b;
    mask_a.SetBit(1);
    mask_b.SetBit(100);

    mask_a.Merge(mask_b);
    CHECK(mask_a.TestBit(1));
    CHECK(mask_a.TestBit(100));
    CHECK(mask_a.PopCount() == 2);
  }

  SECTION("Merge into empty mask copies source") {
    SampleMask empty;
    SampleMask source;
    source.SetBit(0);
    source.SetBit(42);

    empty.Merge(source);
    CHECK(empty.TestBit(0));
    CHECK(empty.TestBit(42));
    CHECK(empty.PopCount() == 2);
  }
}

TEST_CASE("SampleMask::PopCount", "[lancet][cbdg][SampleMask]") {
  SampleMask mask;

  SECTION("Empty mask has zero pop count") {
    CHECK(mask.PopCount() == 0);
  }

  SECTION("PopCount with 10 bits set") {
    for (usize bit_idx = 0; bit_idx < 10; ++bit_idx) {
      mask.SetBit(bit_idx);
    }
    CHECK(mask.PopCount() == 10);
  }

  SECTION("PopCount spanning multiple words") {
    for (usize bit_idx = 0; bit_idx < 128; ++bit_idx) {
      mask.SetBit(bit_idx);
    }
    CHECK(mask.PopCount() == 128);
  }
}

TEST_CASE("SampleMask equality operators", "[lancet][cbdg][SampleMask]") {
  SampleMask mask_a;
  SampleMask mask_b;

  SECTION("Two empty masks are equal") {
    CHECK(mask_a == mask_b);
  }

  SECTION("Same bits set are equal") {
    mask_a.SetBit(5);
    mask_a.SetBit(10);
    mask_b.SetBit(5);
    mask_b.SetBit(10);
    CHECK(mask_a == mask_b);
  }

  SECTION("Different bits set are not equal") {
    mask_a.SetBit(5);
    mask_b.SetBit(6);
    CHECK(mask_a != mask_b);
  }
}

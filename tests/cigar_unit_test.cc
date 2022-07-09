#include "catch2/catch_test_macros.hpp"
#include "lancet2/cigar.h"

TEST_CASE("Can build and operate CigarUnit with CigarOps", "[lancet2::CigarUnit]") {
  SECTION("M CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::ALIGNMENT_MATCH, 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10M");
    CHECK(unit.Operation == lancet2::CigarOp::ALIGNMENT_MATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("I CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::INSERTION, 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10I");
    CHECK(unit.Operation == lancet2::CigarOp::INSERTION);
    CHECK(unit.Length == 10);
  }

  SECTION("D CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::DELETION, 10);
    CHECK(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10D");
    CHECK(unit.Operation == lancet2::CigarOp::DELETION);
    CHECK(unit.Length == 10);
  }

  SECTION("N CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::REFERENCE_SKIP, 10);
    CHECK(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10N");
    CHECK(unit.Operation == lancet2::CigarOp::REFERENCE_SKIP);
    CHECK(unit.Length == 10);
  }

  SECTION("S CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::SOFT_CLIP, 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10S");
    CHECK(unit.Operation == lancet2::CigarOp::SOFT_CLIP);
    CHECK(unit.Length == 10);
  }

  SECTION("H CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::HARD_CLIP, 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10H");
    CHECK(unit.Operation == lancet2::CigarOp::HARD_CLIP);
    CHECK(unit.Length == 10);
  }

  SECTION("P CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::ALIGNMENT_PAD, 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10P");
    CHECK(unit.Operation == lancet2::CigarOp::ALIGNMENT_PAD);
    CHECK(unit.Length == 10);
  }

  SECTION("= CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::SEQUENCE_MATCH, 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10=");
    CHECK(unit.Operation == lancet2::CigarOp::SEQUENCE_MATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("X CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::SEQUENCE_MISMATCH, 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10X");
    CHECK(unit.Operation == lancet2::CigarOp::SEQUENCE_MISMATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("? CigarOp") {
    lancet2::CigarUnit unit(lancet2::CigarOp::UNKNOWN_OP, 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10?");
    CHECK(unit.Operation == lancet2::CigarOp::UNKNOWN_OP);
    CHECK(unit.Length == 10);
  }
}

TEST_CASE("Can build and operate CigarUnit with character ops", "[lancet2::CigarUnit]") {
  SECTION("M CigarOp") {
    lancet2::CigarUnit unit('M', 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10M");
    CHECK(unit.Operation == lancet2::CigarOp::ALIGNMENT_MATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("I CigarOp") {
    lancet2::CigarUnit unit('I', 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10I");
    CHECK(unit.Operation == lancet2::CigarOp::INSERTION);
    CHECK(unit.Length == 10);
  }

  SECTION("D CigarOp") {
    lancet2::CigarUnit unit('D', 10);
    CHECK(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10D");
    CHECK(unit.Operation == lancet2::CigarOp::DELETION);
    CHECK(unit.Length == 10);
  }

  SECTION("N CigarOp") {
    lancet2::CigarUnit unit('N', 10);
    CHECK(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10N");
    CHECK(unit.Operation == lancet2::CigarOp::REFERENCE_SKIP);
    CHECK(unit.Length == 10);
  }

  SECTION("S CigarOp") {
    lancet2::CigarUnit unit('S', 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10S");
    CHECK(unit.Operation == lancet2::CigarOp::SOFT_CLIP);
    CHECK(unit.Length == 10);
  }

  SECTION("H CigarOp") {
    lancet2::CigarUnit unit('H', 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10H");
    CHECK(unit.Operation == lancet2::CigarOp::HARD_CLIP);
    CHECK(unit.Length == 10);
  }

  SECTION("P CigarOp") {
    lancet2::CigarUnit unit('P', 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10P");
    CHECK(unit.Operation == lancet2::CigarOp::ALIGNMENT_PAD);
    CHECK(unit.Length == 10);
  }

  SECTION("= CigarOp") {
    lancet2::CigarUnit unit('=', 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10=");
    CHECK(unit.Operation == lancet2::CigarOp::SEQUENCE_MATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("X CigarOp") {
    lancet2::CigarUnit unit('X', 10);
    CHECK(unit.ConsumesReference());
    CHECK(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10X");
    CHECK(unit.Operation == lancet2::CigarOp::SEQUENCE_MISMATCH);
    CHECK(unit.Length == 10);
  }

  SECTION("? CigarOp") {
    lancet2::CigarUnit unit('a', 10);
    CHECK_FALSE(unit.ConsumesReference());
    CHECK_FALSE(unit.ConsumesQuery());
    CHECK(unit.ToString() == "10?");
    CHECK(unit.Operation == lancet2::CigarOp::UNKNOWN_OP);
    CHECK(unit.Length == 10);
  }
}

TEST_CASE("Can convert AlignmentCigar to string", "[lancet2::AlignmentCigar]") {
  lancet2::AlignmentCigar cigar;
  cigar.emplace_back(lancet2::CigarUnit('M', 100));
  cigar.emplace_back(lancet2::CigarUnit('S', 23));
  cigar.emplace_back(lancet2::CigarUnit('I', 20));
  CHECK(lancet2::ToString(cigar) == "100M23S20I");
}

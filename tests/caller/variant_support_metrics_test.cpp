#include "lancet/base/types.h"
#include "lancet/caller/per_allele_data.h"
#include "lancet/caller/sample_format_data.h"
#include "lancet/caller/variant_support.h"

#include "catch_amalgamated.hpp"

#include <optional>

namespace lancet::caller::tests {

namespace {
// ============================================================================
// Helper: build a ReadEvidence with only the fields relevant to a specific
// metric under test. All other fields use safe defaults.
// ============================================================================
auto MakeEvidence(AlleleIndex allele, i64 aln_start, u32 own_hap_nm, u32 haplotype_id,
                  u32 rname_hash) -> VariantSupport::ReadEvidence {
  return {
      .mInsertSize = 300,
      .mAlignmentStart = aln_start,
      .mAlnScore = 100.0,
      .mFoldedReadPos = 0.25,
      .mRnameHash = rname_hash,
      .mRefNm = 0,
      .mOwnHapNm = own_hap_nm,
      .mAssignedHaplotypeId = haplotype_id,
      .mAllele = allele,
      .mStrand = Strand::FWD,
      .mBaseQual = 30,
      .mMapQual = 60,
      .mIsSoftClipped = false,
      .mIsProperPair = true,
  };
}
}  // namespace

// ============================================================================
// FSSE Tests
// ============================================================================

TEST_CASE("FSSE returns nullopt with fewer than 3 ALT reads",
          "[lancet][caller][VariantSupport][FSSE]") {
  VariantSupport support;
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1003, 0, 1, 101));
  REQUIRE_FALSE(support.ComputeFSSE().has_value());
}

TEST_CASE("FSSE returns zero when all ALT starts land in the same 3bp bin",
          "[lancet][caller][VariantSupport][FSSE]") {
  VariantSupport support;
  // Positions 999, 1000, 1001 all map to bin 333 (floor(x/3))
  support.AddEvidence(MakeEvidence(1, 999, 0, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 101));
  support.AddEvidence(MakeEvidence(1, 1001, 0, 1, 102));
  support.AddEvidence(MakeEvidence(1, 999, 0, 1, 103));
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 104));

  auto const fsse = support.ComputeFSSE();
  REQUIRE(fsse.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(fsse.value(), Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("FSSE returns high entropy for diverse start positions",
          "[lancet][caller][VariantSupport][FSSE]") {
  VariantSupport support;
  // 10 reads at 100bp spacing — each in a unique 3bp bin
  for (u32 i = 0; i < 10; ++i) {
    support.AddEvidence(MakeEvidence(1, static_cast<i64>(i) * 100, 0, 1, 200 + i));
  }

  auto const fsse = support.ComputeFSSE();
  REQUIRE(fsse.has_value());
  // 10 unique bins, normalized by log2(min(10,20)) = log2(10) ≈ 3.32
  // H = log2(10) ≈ 3.32, normalized = 1.0
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(fsse.value(), Catch::Matchers::WithinAbs(1.0, 1e-6));
}

TEST_CASE("FSSE absorbs exonuclease fraying within 3bp bins",
          "[lancet][caller][VariantSupport][FSSE]") {
  VariantSupport support;
  // 6 reads at 1000, 1001, 1002, 1000, 1001, 1002
  // floor(1000/3)=333, floor(1001/3)=333, floor(1002/3)=334
  // Two bins: 333 (4 reads) and 334 (2 reads)
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1001, 0, 1, 101));
  support.AddEvidence(MakeEvidence(1, 1002, 0, 1, 102));
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 103));
  support.AddEvidence(MakeEvidence(1, 1001, 0, 1, 104));
  support.AddEvidence(MakeEvidence(1, 1002, 0, 1, 105));

  auto const fsse = support.ComputeFSSE();
  REQUIRE(fsse.has_value());
  // Two bins: p1=4/6, p2=2/6. H = -(4/6*log2(4/6) + 2/6*log2(2/6)) ≈ 0.918
  // Normalized by log2(min(6, 20)) = log2(6) ≈ 2.585
  // Normalized H ≈ 0.918 / 2.585 ≈ 0.355
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE(fsse.value() > 0.3);
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE(fsse.value() < 0.4);
}

TEST_CASE("FSSE ignores REF reads — only ALT allele starts contribute",
          "[lancet][caller][VariantSupport][FSSE]") {
  VariantSupport support;
  // 10 REF reads at diverse positions — should be ignored
  for (u32 i = 0; i < 10; ++i) {
    support.AddEvidence(MakeEvidence(0, static_cast<i64>(i) * 100, 0, 0, 300 + i));
  }
  // Only 2 ALT reads — below threshold
  support.AddEvidence(MakeEvidence(1, 500, 0, 1, 400));
  support.AddEvidence(MakeEvidence(1, 600, 0, 1, 401));
  REQUIRE_FALSE(support.ComputeFSSE().has_value());
}

// ============================================================================
// AHDD Tests
// ============================================================================

TEST_CASE("AHDD returns nullopt when REF group is empty",
          "[lancet][caller][VariantSupport][AHDD]") {
  VariantSupport support;
  support.AddEvidence(MakeEvidence(1, 1000, 2, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1003, 3, 1, 101));
  REQUIRE_FALSE(support.ComputeAHDD().has_value());
}

TEST_CASE("AHDD returns nullopt when ALT group is empty",
          "[lancet][caller][VariantSupport][AHDD]") {
  VariantSupport support;
  support.AddEvidence(MakeEvidence(0, 1000, 1, 0, 100));
  support.AddEvidence(MakeEvidence(0, 1003, 0, 0, 101));
  REQUIRE_FALSE(support.ComputeAHDD().has_value());
}

TEST_CASE("AHDD returns zero when both groups have equal mean NM",
          "[lancet][caller][VariantSupport][AHDD]") {
  VariantSupport support;
  // REF reads with NM=2 each against their assigned haplotype (REF)
  support.AddEvidence(MakeEvidence(0, 1000, 2, 0, 100));
  support.AddEvidence(MakeEvidence(0, 1003, 2, 0, 101));
  // ALT reads with NM=2 each against their assigned haplotype (ALT)
  support.AddEvidence(MakeEvidence(1, 1006, 2, 1, 102));
  support.AddEvidence(MakeEvidence(1, 1009, 2, 1, 103));

  auto const ahdd = support.ComputeAHDD();
  REQUIRE(ahdd.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(ahdd.value(), Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("AHDD is positive when ALT reads fit their haplotype worse than REF",
          "[lancet][caller][VariantSupport][AHDD]") {
  VariantSupport support;
  // REF reads: NM=1 against REF haplotype (good fit)
  support.AddEvidence(MakeEvidence(0, 1000, 1, 0, 100));
  support.AddEvidence(MakeEvidence(0, 1003, 1, 0, 101));
  // ALT reads: NM=5 against ALT haplotype (poor fit — assembly hallucination)
  support.AddEvidence(MakeEvidence(1, 1006, 5, 1, 102));
  support.AddEvidence(MakeEvidence(1, 1009, 5, 1, 103));

  auto const ahdd = support.ComputeAHDD();
  REQUIRE(ahdd.has_value());
  // AHDD = mean(5,5) - mean(1,1) = 5.0 - 1.0 = 4.0
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(ahdd.value(), Catch::Matchers::WithinAbs(4.0, 1e-6));
}

// ============================================================================
// HSE Tests
// ============================================================================

TEST_CASE("HSE returns nullopt with single haplotype", "[lancet][caller][VariantSupport][HSE]") {
  VariantSupport support;
  for (u32 i = 0; i < 5; ++i) {
    support.AddEvidence(MakeEvidence(1, 1000 + static_cast<i64>(i), 0, 1, 100 + i));
  }
  REQUIRE_FALSE(support.ComputeHSE(1).has_value());
}

TEST_CASE("HSE returns nullopt with fewer than 3 ALT reads",
          "[lancet][caller][VariantSupport][HSE]") {
  VariantSupport support;
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1003, 0, 2, 101));
  REQUIRE_FALSE(support.ComputeHSE(3).has_value());
}

TEST_CASE("HSE returns zero when all ALT reads assigned to same path",
          "[lancet][caller][VariantSupport][HSE]") {
  VariantSupport support;
  // 5 ALT reads, all on haplotype 1
  for (u32 i = 0; i < 5; ++i) {
    support.AddEvidence(MakeEvidence(1, 1000 + static_cast<i64>(i), 0, 1, 100 + i));
  }

  auto const hse = support.ComputeHSE(3);
  REQUIRE(hse.has_value());
  // Only 1 bin → H = 0.0
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(hse.value(), Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("HSE returns 1.0 when ALT reads uniformly split across all haplotypes",
          "[lancet][caller][VariantSupport][HSE]") {
  VariantSupport support;
  // 3 ALT reads, each on a different haplotype (1, 2, 3)
  support.AddEvidence(MakeEvidence(1, 1000, 0, 1, 100));
  support.AddEvidence(MakeEvidence(1, 1003, 0, 2, 101));
  support.AddEvidence(MakeEvidence(1, 1006, 0, 3, 102));

  auto const hse = support.ComputeHSE(3);
  REQUIRE(hse.has_value());
  // H = log2(3), normalized by log2(3) = 1.0
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(hse.value(), Catch::Matchers::WithinAbs(1.0, 1e-6));
}

TEST_CASE("HSE ignores REF reads — only ALT haplotype IDs contribute",
          "[lancet][caller][VariantSupport][HSE]") {
  VariantSupport support;
  // REF reads on various haplotypes — should be ignored
  support.AddEvidence(MakeEvidence(0, 1000, 0, 0, 100));
  support.AddEvidence(MakeEvidence(0, 1003, 0, 1, 101));
  support.AddEvidence(MakeEvidence(0, 1006, 0, 2, 102));
  // Only 2 ALT reads — below threshold
  support.AddEvidence(MakeEvidence(1, 1009, 0, 1, 103));
  support.AddEvidence(MakeEvidence(1, 1012, 0, 2, 104));
  REQUIRE_FALSE(support.ComputeHSE(3).has_value());
}

// ============================================================================
// PDCV — integration testing happens via variant_call_test.cpp since PDCV
// is computed from Graph::Path metadata, not from VariantSupport directly.
// These tests verify the SampleFormatData rendering path.
// ============================================================================

TEST_CASE("PDCV returns nullopt when absent in SampleFormatData",
          "[lancet][caller][VariantCall][PDCV]") {
  SampleFormatData const sample;
  // PDCV defaults to absent (no SetField call)
  REQUIRE_FALSE(sample.GetField(SampleFormatData::PATH_DEPTH_CV).has_value());
}

TEST_CASE("PDCV stores f32 value when present", "[lancet][caller][VariantCall][PDCV]") {
  SampleFormatData sample;
  sample.SetField(SampleFormatData::PATH_DEPTH_CV, 0.45);
  auto const pdcv = sample.GetField(SampleFormatData::PATH_DEPTH_CV);
  REQUIRE(pdcv.has_value());
  // engaged-optional asserted on the prior REQUIRE; clang-tidy does not see through Catch2 macros.
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  REQUIRE_THAT(pdcv.value(), Catch::Matchers::WithinAbs(0.45F, 1e-4F));
}

}  // namespace lancet::caller::tests

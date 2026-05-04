---
description: Lancet2 base/ layer rules — numerical stability, SIMD primitives, statistical effect-size tests, project-wide type aliases. Load when editing src/lancet/base/.
paths:
  - "src/lancet/base/**"
---

# base/ layer rules

`base/` carries the project-wide primitives that every other layer
depends on. Bugs here propagate everywhere; performance regressions
here affect every window in every run.

## Type aliases are project-wide and load-bearing

`base/types.h` defines `i8/i16/i32/i64`, `u8/u16/u32/u64`, `usize`,
`f32`, `f64`. All Lancet2 code uses these aliases instead of
`std::int32_t`/etc. directly. **These aliases are exempt from the
project's snake_case identifier rule** — `validate_cpp_identifiers.py` carries
this exemption explicitly. New aliases added to `types.h` should be
short and follow the same pattern.

## Numerical stability is the base layer's responsibility

`compute_stats.h::OnlineStats` uses Welford's recurrence (m1, m2)
specifically because `Var = E[x²] − (E[x])²` suffers catastrophic
cancellation when the mean is large relative to the spread (e.g., base
qualities mean ≈ 35, σ ≈ 2 — naive formula loses 4+ digits of
precision). Welford's tracks deviations from the running mean so the
subtraction is between similar-magnitude quantities. Never replace
`OnlineStats::Add` with a naive Σx²/n computation.

`OnlineStats::Merge` implements Chan et al.'s parallel merge formula:
`m2_combined = m2₁ + m2₂ + δ²·n₁·n₂/(n₁+n₂)` where `δ = m1₂ − m1₁`.
Per-thread accumulators merge correctly without re-reading data; this
is required for the multi-threaded coverage stats path. A simple
`m2 += other.m2` is silently wrong.

## Mann-Whitney effect size: Z/√N, not raw Z

`mann_whitney.h::MannWhitneyEffectSize` returns Z/√N, not raw Z. The
raw test statistic Z scales with √N by the CLT — at 2000× coverage,
even a trivial 0.5-unit quality difference produces a "highly
significant" Z. Dividing by √N cancels this inflation: the same
biological bias produces the same effect size at 20× and 2000×. ML
features computed from this must be coverage-invariant.

The `std::optional` return type encodes a meaningful distinction:
`std::nullopt` = test cannot run (one group empty — homozygous,
allelic dropout). `0.0` = test ran, found no bias (genuine zero).
Downstream, `std::nullopt` becomes `.` in VCF (untestable); `0.0` is
emitted as `0.0`. Conflating these two is a schema bug — see
`vcf-validator` for the missing-value convention.

Tie correction uses Lehmann's variance formula: `Var(U) = (m·n/12) ·
[(N+1) − Σₖ(tₖ³−tₖ)/(N(N−1))]`. Drop the tie correction term and
discrete data (integer MAPQs, integer base qualities) gets a wrong
variance, biasing the effect size.

## Polar coordinates encode independence, not just transformation

`polar_coords.h` exists because DP and VAF are NOT independent ML
features. DP = AD_Ref + AD_Alt couples the two counts (somatic
AD=(95,5) and germline AD=(50,50) both give DP=100); VAF has
depth-dependent precision. The polar transform produces:

- `PANG = atan2(AD_Alt, AD_Ref)` — identity, depth-independent angle
- `PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²))` — magnitude, log-compressed

These two axes are statistically independent. Any change here that
accidentally couples them (e.g., normalizing PRAD by PANG) defeats
the purpose. The features must remain orthogonal.

## SIMD is platform-conditional, with verified tail handling

`repeat.cpp::IsWithinHammingDist` has three paths: AVX2 on x86
(`_mm256_cmpeq_epi8` + `movemask` + `POPCNT`), NEON on ARM64
(`vceqq_u8` + `vbicq_u8` + `vaddvq_u8` — one cycle vs. invert+sum),
and auto-vectorized fallback. CMake guarantees `-march=x86-64-v3`
(AVX2/BMI2/POPCNT baseline) on x86 and aarch64 baseline NEON on ARM,
so the architecture macros are reliable. **Do not** add a scalar tail
loop after the SIMD body — the existing code uses overlapping
unaligned loads for the tail, which avoids branchy cleanup. New SIMD
additions should follow the same pattern.

The early-exit in `IsWithinHammingDist` after each SIMD chunk is what
collapses the per-pair cost in the O(n²) repeat-detection loop from
O(L) to effectively O(1) for random DNA. Removing the early-exit
silently makes repeat detection 100× slower without breaking
correctness.

## Longdust scoring: GC-bias correction must be GLOBAL, not local

`longdust_scorer.h` implements Q(x) complexity from Li 2025. The
`gc_frac` parameter MUST be a global (genome-wide) or broad regional
fraction, never the local GC of the scored window. A poly-A insertion
locally has 0% GC; passing 0.0 makes the scorer expect all-A k-mers
and produce Q(x)≈0, blinding it to the repeat. The default 0.41 is
human genome-wide GC; non-human callers must supply their organism's
genome-wide value.

GC correction only changes `f(ℓ)` precomputation in the constructor;
the hot-path `Score()` method has zero per-call overhead from the
correction.

## Logging is async with mandatory backpressure

`logging.h::RegisterLancetLogger` initializes spdlog with
`async_overflow_policy::block` and `QUEUE_SIZE = 32768`. Block (not
discard) is intentional: dropped log lines during a 40-hour pipeline
run make post-mortem debugging impossible. Use `LOG_TRACE`/`DEBUG`/
`INFO`/`WARN`/`ERROR`/`CRITICAL` macros — never `spdlog::info` etc.
directly. `LOG_TRACE` is compiled out unless `LANCET_VERBOSE_MODE` is
defined.

## LANCET_ASSERT is debug-only and throws

`assert.h::LANCET_ASSERT(condition)` expands to a throw with
`std::source_location` only when `LANCET_DEBUG_MODE` is defined; in
release it expands to `((void)0)`. Conditions in assertions must be
side-effect-free. Asserts that do real work in debug but compile away
in release are the worst-class bug — treat the macro like
`assert(3)`.

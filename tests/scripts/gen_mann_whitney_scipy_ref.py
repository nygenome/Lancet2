"""Regenerate the Mann-Whitney effect-size reference TSV using scipy.

This is a one-shot generator for `tests/data/base/mann_whitney_scipy_ref.tsv`,
the golden-master fixture that
`tests/base/mann_whitney_test.cpp` cross-validates the project's
`MannWhitneyEffectSize` against. Per the project's test-fixture-script
convention, all generators live under `tests/scripts/` while their
generated outputs live under `tests/data/<layer>/`.

Project's MannWhitneyEffectSize semantics
-----------------------------------------
The in-tree implementation in `src/lancet/base/mann_whitney.h` returns the
coverage-normalized effect size `Z / sqrt(N)`, where:

  - U      = sum of mid-ranks for the alt sample minus n_alt*(n_alt+1)/2
  - E[U]   = n_ref * n_alt / 2 under H0
  - Var(U) = (n_ref*n_alt/12) * [(N+1) - sum(t^3 - t) / (N*(N-1))]
             with Lehmann's tie correction
  - Z      = (U - E[U]) / sqrt(Var(U))     (no continuity correction)
  - eff    = Z / sqrt(N)                   where N = n_ref + n_alt
  - Empty group     -> nullopt (encoded as `nan` in the TSV; the C++ test
                                  treats `nan` as "test untestable").
  - Zero variance   -> 0.0     (genuine zero, returned when all values
                                  are identical or an extreme tie-correction
                                  drives Var(U) to 0).

This script reuses scipy's mannwhitneyu to obtain U (matching the project's
sign convention by treating the alt sample as the first argument), then
computes Var(U) with Lehmann's tie-corrected formula and the effect size
`Z / sqrt(N)` exactly the way `MannWhitneyEffectSize` does.

Regeneration
------------
    pixi run -e test-fixtures \\
        python tests/scripts/gen_mann_whitney_scipy_ref.py

The output `tests/data/base/mann_whitney_scipy_ref.tsv` is committed
alongside this script so test runs do not depend on scipy at runtime.

TSV schema
----------
Tab-separated, one row per case, with a single header line. Columns:
    n_ref               int    Size of the REF sample.
    n_alt               int    Size of the ALT sample.
    groups_seed         int    Seed used by numpy's Generator to produce the
                                two samples. Documented for reproducibility;
                                C++ reads ref_vals/alt_vals directly rather
                                than re-seeding (numpy's PCG64 is not
                                available in C++).
    ref_vals            str    Comma-separated f64 values for the REF sample.
                                Empty string when n_ref == 0.
    alt_vals            str    Comma-separated f64 values for the ALT sample.
                                Empty string when n_alt == 0.
    expected_effect_size  f64  Project's `Z / sqrt(N)` effect size, or `nan`
                                when n_ref == 0 or n_alt == 0 (untestable).

The C++ test uses `Catch::Approx().margin(...)` to absorb f64 rounding
between scipy's path and the in-tree path.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

# This script lives in tests/scripts/. The output TSV lives next to its
# consumer in tests/data/base/. Resolve the project root by going two levels
# up from this script's location, then descending to the fixture path.
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_TSV = PROJECT_ROOT / "tests" / "data" / "base" / "mann_whitney_scipy_ref.tsv"


def project_effect_size(ref: np.ndarray, alt: np.ndarray) -> float:
    """Reproduce src/lancet/base/mann_whitney.h::MannWhitneyEffectSize.

    Empty groups -> NaN. Zero variance -> 0.0. Otherwise: Z/sqrt(N) using
    Lehmann's tie-corrected variance and no continuity correction.
    """
    if len(ref) == 0 or len(alt) == 0:
        return float("nan")

    n_ref = float(len(ref))
    n_alt = float(len(alt))
    n_total = n_ref + n_alt

    # scipy returns U for the FIRST argument. Pass alt first to match the
    # project's U_alt sign convention (U_alt = R_alt - n_alt*(n_alt+1)/2,
    # which is what scipy.mannwhitneyu(alt, ref).statistic returns when
    # use_continuity=False).
    u_alt = mannwhitneyu(alt, ref, alternative="two-sided",
                         use_continuity=False, method="asymptotic").statistic

    mean_u = n_ref * n_alt / 2.0

    # Tie correction: sum of (t^3 - t) over each tie group of size t.
    pooled = np.concatenate([ref, alt])
    counts = Counter(pooled.tolist())
    tie_correction = sum((t * t * t) - t for t in counts.values())

    var_u = (n_ref * n_alt / 12.0) * (
        (n_total + 1.0) - (tie_correction / (n_total * (n_total - 1.0)))
    )
    if var_u <= 0.0:
        return 0.0

    z_score = (u_alt - mean_u) / math.sqrt(var_u)
    return z_score / math.sqrt(n_total)


def synthesize(seed: int, n_ref: int, n_alt: int, alt_shift: float
               ) -> tuple[np.ndarray, np.ndarray]:
    """Two integer-valued samples drawn from a 60-base-quality-like
    distribution with an optional shift on the ALT side.

    Integer values exercise the rank-tie path; the shift produces a
    detectable effect at most sample sizes."""
    rng = np.random.default_rng(seed)
    ref = rng.integers(low=20, high=60, size=n_ref).astype(np.float64)
    alt = (rng.integers(low=20, high=60, size=n_alt).astype(np.float64)
           + alt_shift)
    return ref, alt


def fmt_f64(value: object) -> str:
    """Round-trip-safe f64 string. `repr(np.float64(x))` in numpy >=2.0
    wraps the value (\"np.float64(...)\"); `repr(float(x))` does not, and
    yields the shortest decimal that round-trips through std::stod."""
    return repr(float(value))


def write_tsv(rows: Iterable[dict]) -> None:
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ("n_ref", "n_alt", "groups_seed", "ref_vals", "alt_vals",
            "expected_effect_size")
    with OUTPUT_TSV.open("w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            ref_str = ",".join(fmt_f64(v) for v in row["ref_vals"])
            alt_str = ",".join(fmt_f64(v) for v in row["alt_vals"])
            out.write(
                f"{row['n_ref']}\t{row['n_alt']}\t{row['groups_seed']}\t"
                f"{ref_str}\t{alt_str}\t{fmt_f64(row['expected_effect_size'])}\n"
            )


def build_corpus() -> list[dict]:
    """Hand-curated cases covering the input dimensions the project's
    MannWhitneyEffectSize must handle correctly:

      - Balanced small / medium sizes
      - Asymmetric sample sizes (1 vs n, 9 vs 1)
      - With and without an ALT shift (shift=0 -> near-zero effect)
      - With heavy ties (low-cardinality integer values)
      - With a single-element group on either side
    """
    cases: list[dict] = []

    config: list[tuple[int, int, int, float]] = [
        (5, 5, 1001, 0.0),       # balanced, no shift
        (5, 5, 1001, 5.0),       # balanced, mild shift
        (10, 10, 1002, 0.0),     # balanced, no shift
        (10, 10, 1002, 10.0),    # balanced, large shift
        (20, 20, 1003, 5.0),     # balanced, mild shift
        (50, 50, 1004, 2.0),     # bigger sample, small shift
        (100, 100, 1005, 1.0),   # large balanced, faint shift
        (10, 1, 1006, 5.0),      # 10 ref vs 1 alt, shifted
        (1, 10, 1007, -5.0),     # 1 ref vs 10 alt, shifted negative
        (9, 1, 1008, 0.0),       # asymmetric, no shift
        (5, 30, 1009, 3.0),      # heavily asymmetric
        (3, 3, 1010, 0.0),       # tiny balanced, no shift
        (3, 3, 1010, 4.0),       # tiny balanced, shift
        (15, 15, 1011, 0.5),     # small shift
        (50, 5, 1012, 5.0),      # asymmetric large/small with shift
    ]

    for n_ref, n_alt, seed, shift in config:
        ref, alt = synthesize(seed, n_ref, n_alt, shift)
        cases.append({
            "n_ref": n_ref,
            "n_alt": n_alt,
            "groups_seed": seed,
            "ref_vals": ref,
            "alt_vals": alt,
            "expected_effect_size": project_effect_size(ref, alt),
        })

    # Empty-group cases (untestable; effect size is NaN).
    cases.append({"n_ref": 0, "n_alt": 5, "groups_seed": 0,
                  "ref_vals": np.array([], dtype=np.float64),
                  "alt_vals": np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
                  "expected_effect_size": float("nan")})
    cases.append({"n_ref": 5, "n_alt": 0, "groups_seed": 0,
                  "ref_vals": np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
                  "alt_vals": np.array([], dtype=np.float64),
                  "expected_effect_size": float("nan")})

    # All-tied case (effect size is genuinely 0).
    tied = np.full(5, 7.0)
    cases.append({"n_ref": 5, "n_alt": 5, "groups_seed": 0,
                  "ref_vals": tied, "alt_vals": tied,
                  "expected_effect_size": 0.0})

    return cases


def main() -> None:
    rows = build_corpus()
    write_tsv(rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT_TSV}")


if __name__ == "__main__":
    main()

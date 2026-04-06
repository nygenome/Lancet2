# Polar Coordinate Features (PRAD / PANG)

These per-sample FORMAT fields transform Cartesian allele depth coordinates
(`AD_Ref`, `AD_Alt`) into polar coordinates, orthogonalizing biological
identity from sequencing depth for downstream ML variant classification.

---

## Motivation

When plotting `AD_Ref` (X) vs `AD_Alt` (Y), variant classes suffer two
geometric distortions in Cartesian space:

1. **The "Diagonal" Problem**: Biologically identical variants (e.g.,
   heterozygous germline at 50% VAF) form elongated diagonals spanning
   `(5,5)` at low depth to `(500,500)` at high depth. Distance-based ML
   algorithms (K-Means, KNN) misinterpret opposite ends of this diagonal
   as different classes.

2. **The "Wedge" (Fan) Problem**: Variance in allele counts scales with
   depth. A 50% VAF at 10× has tight variance; at 1000× it fans out.
   This creates wedge-shaped clusters whose width constantly changes,
   preventing clean decision boundaries.

The polar transform "unrolls" the wedge and "collapses" the diagonal into
a uniform rectangular strip where ML models can draw simple flat boundaries.

---

## Polar Radius (`PRAD`)

**Computation**: `PRAD = log10(1 + sqrt(AD_Ref² + AD_Alt²))`, computed per-sample.

The log10-compressed Euclidean magnitude of the allele depth vector. Encodes
total signal strength (sequencing depth) orthogonal to allele fraction, with
built-in normalization for cross-coverage ML generalization.

**Why log10?** The raw Euclidean radius `sqrt(AD_Ref² + AD_Alt²)` scales
linearly with coverage, spanning 0 to >2800 at 2000× depth — a 100× dynamic
range. An ML model trained at 30× WGS would see completely different PRAD
distributions when applied to 100× WES data. The `log10(1+r)` compression
reduces this to a compact [0, ~3.5] range while preserving monotonicity:

| Coverage | Het 50% VAF | PRAD (log10) | Interpretation |
|:---------|:------------|:-------------|:---------------|
| 0× | — | 0.0 | No reads — no evidence |
| 20× | AD=(10,10) | 1.15 | Low confidence |
| 60× | AD=(30,30) | 1.63 | Moderate confidence |
| 100× | AD=(50,50) | 1.85 | Good confidence |
| 1000× | AD=(500,500) | 2.85 | Very high confidence |
| 2000× | AD=(1000,1000) | 3.15 | Extreme confidence |

**Value range**: [0, ~3.5]

**Interpretation**: Higher PRAD means more sequencing evidence, regardless
of whether it supports REF or ALT. Acts as a confidence tie-breaker for
low-angle variants:

| PRAD | Meaning |
|:-----|:--------|
| < 1.0 + low PANG | Stochastic noise — insufficient evidence |
| < 1.0 + high PANG | Low-coverage heterozygous/homozygous call |
| > 2.0 + low PANG | True low-frequency variant (somatic/mosaic) with deep sequencing support |
| > 2.0 + high PANG | Well-supported heterozygous/homozygous variant |

**Coverage stability**: PRAD still varies with coverage (1.15 at 20× to 3.15
at 2000×), but the log10 compression reduces the dynamic range from 100× to
3× and preserves monotonic ordering. For ML models that use PRAD as a
confidence weight, this bounded range avoids feature-scale domination that
occurs with raw allele counts. Tree-based models (XGBoost, Random Forest)
are invariant to monotonic transforms and handle this naturally; linear
models benefit from the compressed range.

---

## Polar Angle (`PANG`)

**Computation**: `PANG = atan2(AD_Alt, AD_Ref)`, computed per-sample.

The arctangent of the allele depth ratio, in radians. Encodes the biological
identity of the variant — the allele fraction — independent of total depth.

**Value range**: [0, π/2] radians (≈ [0, 1.571])

**Interpretation**:

| PANG (radians) | PANG (degrees) | Meaning |
|:---------------|:---------------|:--------|
| ≈ 0.0 | 0° | Homozygous REF — no ALT evidence |
| < 0.524 | < 30° | Low-frequency signal (somatic, mosaic, or artifact) |
| ≈ 0.785 | 45° | Heterozygous germline — balanced REF/ALT |
| ≈ 1.571 | 90° | Homozygous ALT — nearly all reads are ALT |

**Coverage stability**: Perfectly coverage-invariant. PANG is a pure ratio
(`atan2(AD_Alt, AD_Ref)`) — identical allele fractions produce identical
angles at any depth. A 50% VAF het produces PANG ≈ 0.785 (45°) whether
the sample is 20× or 2000×. This is the key advantage of the polar transform:
allele identity is algebraically separated from sequencing depth.

---

## Tumor-Normal "Four-Feature" Paradigm

Each sample independently computes its own PRAD and PANG from its own allele
depths. The downstream ML model receives a four-element feature vector:

`[PRAD_Normal, PANG_Normal, PRAD_Tumor, PANG_Tumor]`

Both PRAD (log10-compressed, range [0, ~3.5]) and PANG (ratio-based, range
[0, π/2]) are coverage-invariant, ensuring the four-feature paradigm produces
comparable values regardless of whether the sample is 20× or 2000×.

This avoids cross-sample relational metrics that would violate VCF FORMAT
semantics and eliminates sentinel/NaN imputation issues. The model learns
the relational delta itself:

| Pattern | PANG Normal | PANG Tumor | Classification |
|:--------|:------------|:-----------|:---------------|
| Germline | ≈ π/4 (45°) | ≈ π/4 (45°) | Inherited variant — same allele fraction in both |
| Somatic | ≈ 0 (0°) | ≈ π/6..π/4 | Tumor-specific variant — normal has no ALT support |
| LOH | ≈ π/4 (45°) | ≈ π/2 (90°) | Loss of heterozygosity — tumor lost REF allele |

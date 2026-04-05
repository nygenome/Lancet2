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

**Computation**: `PRAD = sqrt(AD_Ref² + AD_Alt²)`, computed per-sample.

The Euclidean magnitude of the allele depth vector. Encodes total signal
strength (sequencing depth) orthogonal to allele fraction.

**Value range**: [0, ∞)

**Interpretation**: Higher PRAD means more sequencing evidence, regardless
of whether it supports REF or ALT. Acts as a confidence tie-breaker for
low-angle variants:

| PRAD | Meaning |
|:-----|:--------|
| Low + low PANG | Stochastic noise — insufficient evidence |
| Low + high PANG | Low-coverage heterozygous/homozygous call |
| High + low PANG | True low-frequency variant (somatic/mosaic) with deep sequencing support |
| High + high PANG | Well-supported heterozygous/homozygous variant |

**Note for ML pipelines**: PRAD scales from 0 to >10,000 depending on read
depth. For neural networks, apply `log10(PRAD + 1)` or `StandardScaler`
before training. Tree-based models (XGBoost, Random Forest) are
scale-invariant and can use PRAD directly.

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

---

## Tumor-Normal "Four-Feature" Paradigm

Each sample independently computes its own PRAD and PANG from its own allele
depths. The downstream ML model receives a four-element feature vector:

`[PRAD_Normal, PANG_Normal, PRAD_Tumor, PANG_Tumor]`

This avoids cross-sample relational metrics that would violate VCF FORMAT
semantics and eliminates sentinel/NaN imputation issues. The model learns
the relational delta itself:

| Pattern | PANG Normal | PANG Tumor | Classification |
|:--------|:------------|:-----------|:---------------|
| Germline | ≈ π/4 (45°) | ≈ π/4 (45°) | Inherited variant — same allele fraction in both |
| Somatic | ≈ 0 (0°) | ≈ π/6..π/4 | Tumor-specific variant — normal has no ALT support |
| LOH | ≈ π/4 (45°) | ≈ π/2 (90°) | Loss of heterozygosity — tumor lost REF allele |

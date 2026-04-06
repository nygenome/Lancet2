# Sequence Complexity (`SEQ_CX`)

This optional INFO field captures 11 orthogonal sequence complexity features
distilled from raw multi-scale metrics (homopolymer runs, Shannon entropy,
LongdustQ k-mer concentration, tandem repeat motifs). It requires
`--enable-sequence-complexity-features` to enable.

All 11 features are computed from assembled haplotype sequences and are
**completely coverage-invariant**: they measure properties of the genomic
sequence itself, not properties of the reads. A poly-A run is 12bp long
regardless of whether 20 or 2000 reads cover it.

**INFO tag**: `SEQ_CX` (11 comma-separated values)

---

## Design: Context vs. Perturbation Paradigm

Raw sequence complexity metrics at multiple scales (±5/10/20/50/100bp, full
haplotype) are highly collinear — HRun at ±5bp, ±10bp, and ±20bp are ~0.95
correlated. Additive ML models (EBMs/GAMs) split importance weight across
correlated features, producing noisy shape functions instead of confident ones.

The distillation separates metrics into three orthogonal groups:

1. **Context** (REF-only): "How brittle is the genome here, regardless of the variant?"
2. **Deltas** (ALT−REF): "How did the variant alter the local complexity?"
3. **TR Motif** (ALT-only): "What is the repeat environment of the final allele?"

This enables the model to learn rescue logic: even if context is dangerous
(high ContextHRun), a negative DeltaHRun (variant broke the homopolymer)
rescues the call.

---

## Feature Layout

### A. Context (4 features — strictly REF)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 0 | ContextHRun | Integer | [0, ∞) | Max homopolymer run in REF ±20bp. A poly-A of length 12 means any 1bp INDEL within 20bp is a stutter candidate. |
| 1 | ContextEntropy | Float | [0.0, 2.0] | Shannon entropy H = −Σ pᵢ log₂(pᵢ) in REF ±20bp. 0 = single base type (homopolymer), 2.0 = perfectly balanced ACGT. Low entropy ≈ low sequence diversity ≈ harder to map uniquely. |
| 2 | ContextFlankLQ | Float | [0.0, ~1.6] | log₁p-squashed LongdustQ (k=4) in REF ±50bp. Captures microsatellite density at a scale where k-mer statistics are reliable. The log₁p transform compresses heavy tails (telomeric LQ > 4.0 → ~1.6). |
| 3 | ContextHaplotypeLQ | Float | [0.0, ~1.6] | log₁p-squashed LongdustQ (k=7) on full REF haplotype. Captures macro-scale repetitive structure of the assembled contig. Only REF because a 5bp INDEL changes < 0.5% of a 1000bp k-mer distribution. |

### B. Deltas (3 features — ALT minus REF)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 4 | DeltaHRun | Integer | [−∞, +∞) | ALT ±5bp HRun − REF ±5bp HRun. Positive = variant extended a homopolymer (artifact signal). Negative = variant broke a homopolymer (rescue signal). Computed at ±5bp to maximize sensitivity to immediate allelic change. |
| 5 | DeltaEntropy | Float | [−2.0, +2.0] | ALT ±10bp entropy − REF ±10bp entropy. Negative = variant reduced sequence diversity (gap mimicking deletion). Computed at ±10bp to balance between ultra-local and micro sensitivity. |
| 6 | DeltaFlankLQ | Float | [−1.6, +1.6] | log₁p(ALT ±50bp LQ) − log₁p(REF ±50bp LQ). Positive = variant exacerbated microsatellite repetitiveness. Computed in log-space to normalize the heavy-tailed LQ distribution. |

### C. TR Motif (4 features — strictly ALT ±50bp)

| Index | Field | Type | Range | Description |
|:------|:------|:-----|:------|:------------|
| 7 | TrAffinity | Float | [0.0, 1.0] | 1/(1+dist) where dist = distance to nearest tandem repeat in ALT ±50bp. Sentinel-safe: dist < 0 (no TR found) → 0.0. dist = 0 (sitting on TR) → 1.0. Monotonically decaying with distance. |
| 8 | TrPurity | Float | [0.0, 1.0] | Purity (1 − errors/span) of nearest TR. dist < 0 → 0.0. A purity of 0.95 means near-perfect repeat fidelity, strongly predicting polymerase stutter. |
| 9 | TrPeriod | Integer | [0, 6] | Period of nearest TR. dist < 0 → 0. Period dictates the biophysics of slippage: period 1 (homopolymer) has highest stutter rate, period 6 (hexanucleotide) has lowest. |
| 10 | IsStutterIndel | Integer | {0, 1} | 1 if INDEL length ≤ repeat period AND variant is adjacent to a TR (dist ≤ 1bp). This is the canonical Illumina polymerase stutter signature. |

---

## Transform Details

### Sentinel-Safe TR Features

Raw motif detection uses `dist_to_nearest_tr = −1` as a sentinel for "no TR
found." This breaks EBM numerical binning because −1 sorts between −2 and 0,
making "no TR" appear closer to "sitting on the repeat" than a distance of 10.

The `TrAffinity` inverse transform (`1/(1+dist)`) maps the sentinel into the
continuous [0, 1] range with 0 = no TR and 1 = overlapping TR. All other TR
features are zeroed when no TR is found.

### Log₁p Squashing

LongdustQ reaches 4.0+ in telomeric and pericentromeric regions, creating
extreme heavy tails. The `log₁p(max(0, x))` transform compresses these tails
while preserving monotonicity and differentiability at zero.

### Multi-Haplotype and Multi-Allelic Merging

When a variant appears on multiple ALT haplotypes, the scorer produces one
`SequenceComplexity` per haplotype and merges via element-wise max (pessimistic
worst-case). Multi-allelic variants at the same locus are similarly merged.

---

## GC-Bias Correction

LongdustQ uses a null model for expected k-mer frequencies. By default, it
assumes GC fraction = 0.41 (human genome-wide average). Without correction,
AT-rich regions are scored as artificially "repetitive."

```bash
--genome-gc-bias 0.41    # Human genome (default)
--genome-gc-bias 0.5     # Uniform model (no correction)
--genome-gc-bias 0.42    # Mouse genome
```

---

## Coverage Stability

All 11 SEQ_CX features are **perfectly coverage-invariant** because they
are computed from assembled haplotype sequences, not from read-level metrics.
The scoring pipeline works as follows:

1. The de Bruijn graph assembler produces haplotype sequences (contigs)
   from the reads in each window.
2. The sequence complexity scorer operates on these haplotype strings.
3. Homopolymer runs, Shannon entropy, LongdustQ, and TR motif detection
   all operate on the sequence characters, not on read counts or qualities.

The only indirect coverage dependency is that the assembler itself requires
sufficient coverage to produce correct haplotypes. Below ~10×, the
assembler may fail to construct the ALT haplotype, in which case no variant
is called and no SEQ_CX is computed. Above ~15×, the assembled haplotypes
are identical regardless of depth, and all 11 features are deterministic:

| Feature | Coverage dependency |
|:--------|:--------------------|
| ContextHRun | None — counts bases in REF string |
| ContextEntropy | None — base frequency ratio in REF string |
| ContextFlankLQ | None — k-mer concentration in REF string |
| ContextHaplotypeLQ | None — k-mer concentration in full REF haplotype |
| DeltaHRun | None — ALT string vs REF string |
| DeltaEntropy | None — ALT string vs REF string |
| DeltaFlankLQ | None — ALT string vs REF string |
| TrAffinity | None — motif distance in ALT string |
| TrPurity | None — motif purity in ALT string |
| TrPeriod | None — motif period in ALT string |
| IsStutterIndel | None — binary flag from ALT string |

---

## Motif Detection Parameters

The tandem repeat motif detector scans for exact and approximate repeats.
Default configuration is optimized for somatic variant calling on Illumina
short reads:

| Parameter | Default | Rationale |
|:----------|:--------|:----------|
| `max_period` | 6 | Covers homopolymers through hexanucleotide repeats. Period 1–3 account for > 90% of Illumina INDEL errors. |
| `min_copies_exact` | 2.5 | Catches ATATAT (3× AT) but not ATAT (2× AT). |
| `min_copies_approx` | 3.0 | Higher threshold for approximate matches to compensate for relaxed matching. |
| `max_edits_per_unit` | 1 | Allows 1 mismatch/indel per repeat unit. |
| `min_purity` | 0.75 | Below 75% purity, the repeat is too degraded to cause polymerase stutter. |

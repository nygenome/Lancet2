# Alignment-Derived Annotations

These per-sample FORMAT fields detect alignment, sequencing, and library
preparation artifacts. All fields are designed to be **coverage-invariant**:
the same biological signal produces the same annotation value regardless of
sequencing depth, enabling ML models trained at one coverage to generalize
across 20×–2000× without retraining.

SCA, FLD, and MQCD use metadata from the original BAM/CRAM alignment.
RPCD, BQCD, and ASMD use metrics computed from the minimap2 re-alignment
during genotyping. SDFC uses window-level BAM coverage.

---

## Soft Clip Asymmetry (`SCA`)

**Purpose**: Detect unresolved larger structural variant events that
masquerade as smaller local variants in the assembly.

**Computation**:

1. During BAM/CRAM ingestion, each read is flagged as "soft-clipped" if
   the total soft-clip bases in its original whole-genome alignment CIGAR
   are ≥ 6% of the read length.
2. After genotyping assigns reads to REF or ALT alleles, the soft-clip
   fraction for each allele is computed:
   - `alt_frac = alt_soft_clip_count / alt_total_count` (0 if no ALT reads)
   - `ref_frac = ref_soft_clip_count / ref_total_count` (0 if no REF reads)
3. `SCA = alt_frac − ref_frac`, computed per-sample.

**Value range**: [−1.0, 1.0]

**Interpretation**:

| SCA Range | Meaning |
|:----------|:--------|
| ≈ 0.0 | Symmetric soft-clipping between ALT and REF — likely no hidden event |
| > 0.1 | ALT reads disproportionately soft-clipped — possible unresolved SV, translocation, or complex event |
| < −0.1 | REF reads more clipped — unusual, may indicate mapping artifact near the reference path |

**Coverage stability**: Inherently coverage-invariant — computed as a ratio
of fractions. A 30% soft-clip rate in ALT reads produces SCA ≈ 0.3 at any
depth. Using SCA together with SDFC helps distinguish real SVs (high SCA,
normal SDFC) from mapping artifacts (high SCA, elevated SDFC).

---

## Fragment Length Delta (`FLD`)

**Purpose**: Detect chimeric library artifacts (artificial bridging) or
somatic cfDNA fragment length shifts.

**Computation**:

1. During BAM/CRAM ingestion, each read's template length (`TLEN` /
   `bam1_t::core.isize`) is captured from the original alignment.
2. After genotyping, for each sample, properly-paired reads with non-zero
   insert size are grouped by their assigned allele.
3. Mean insert sizes are computed for REF-supporting and ALT-supporting
   reads separately.
4. `FLD = |mean_alt_isize − mean_ref_isize|`, computed per-sample.

**Value range**: [0, ∞)

**Interpretation**:

| FLD Range | Meaning |
|:----------|:--------|
| < 20 bp | Normal variation — insert size distributions are consistent |
| 20–100 bp | Moderate discrepancy — inspect manually |
| > 100 bp | Large discrepancy — likely chimeric artifact, library prep issue, or cfDNA fragment size shift |

**Coverage stability**: The mean insert size converges rapidly (by ~20
reads per allele). FLD shows minor variation at very low coverage (N<5 per
allele) due to sampling noise, but is effectively stable above 20×.
Interpret FLD jointly with RPCD: a true structural variant often produces
both elevated FLD and edge-biased RPCD.

---

## Mapping Quality Cohen's D (`MQCD`)

**Purpose**: Detect paralogous mismapping — situations where ALT-supporting
reads originate from a different genomic locus (e.g., a segmental duplication)
and are assigned artificially low mapping quality by the aligner.

**Statistical method**: Coverage-normalized effect size from the Mann-Whitney
U test (Wilcoxon Rank-Sum test). The raw Z-score is divided by √N (where
N = total reads) to remove the √N power amplification from the Central Limit
Theorem, recovering a standardized effect size analogous to Cohen's d.

**Motivation for Z/√N normalization**: The raw Mann-Whitney Z-score scales
with √N: the same mild ALT MAPQ depression produces Z ≈ −1.5 at 20× but
Z ≈ −14.9 at 2000×. This makes raw Z-scores unusable for ML models that
must generalize across coverages. Dividing by √N produces a coverage-invariant
effect size that measures the *magnitude* of the MAPQ difference, not the
statistical significance.

**Computation**:

1. Collect original alignment MAPQ (`bam1_t::core.qual`) for all reads,
   grouped by their genotyper-assigned allele (REF vs ALT).
2. Pool all observations into a single ranked sequence. Tied values
   receive the mean (mid-rank) of the positions they span.
3. Compute the U statistic: `U = R_alt − n_alt·(n_alt+1)/2`.
4. Apply the tie-corrected variance formula (Lehmann, 2006):
   `Var(U) = (m·n/12) · [(N+1) − Σ(tₖ³−tₖ)/(N·(N−1))]`
5. `MQCD = Z / √N = [(U − E[U]) / √Var(U)] / √(m+n)`, computed per-sample.

**Value range**: [−2, +2] typically, or `0.0` (untestable)

**Coverage stability**: A constant mild bias (ALT MAPQ 2 units lower)
produces MQCD ≈ −0.34 at **every** depth from 20× to 2000× (< 3% variation).

**Interpretation**:

| MQCD Range | Meaning |
|:-----------|:--------|
| ≈ 0.0 | No systematic MAPQ difference — strong evidence for true variant. Also returned when the test is untestable (empty REF or ALT group, all identical MAPQs). |
| −0.2 to −0.5 | Moderate ALT MAPQ depression — possible repetitive region, inspect manually |
| < −0.5 | Strong ALT mismapping signal — likely paralogous or multi-mapping artifact |
| > 0 | ALT MAPQ higher than REF — unusual, may indicate REF mapping issues |

**Manual interpretation tip**: Use MQCD together with SDFC. Paralogous
mismapping typically produces both negative MQCD (low ALT MAPQ) and
elevated SDFC (excess depth from collapsed paralogs). A site with
MQCD < −0.3 and SDFC > 2.0 is very likely a false positive from
segmental duplication.

**References**:

- Mann, H.B. & Whitney, D.R. (1947). *Annals of Mathematical Statistics*, 18(1), 50–60.
- Lehmann, E.L. (2006). *Nonparametrics: Statistical Methods Based on Ranks*, Springer.
- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*, 2nd ed.

---

## Read Position Cohen's D (`RPCD`)

**Purpose**: Detect systematic read-edge bias in ALT-supporting reads.
Artifacts from 3' quality degradation and 5' soft-clip misalignment produce
false variant calls that cluster at read extremes.

**Key insight — folded read position**: The test uses the **folded** read
position rather than the raw position. Raw positions are bimodal for
edge-biased artifacts (clustered at both 5' and 3' ends), but the mean
of positions 5 and 145 in a 150bp read ≈ 75 — indistinguishable from a
truly centered variant. Folding maps both ends to the same low-value space:

$$P_{folded} = \min(P_{raw},\; L - 1 - P_{raw})$$

This converts the bimodal trap into a unidirectional signal: "Are ALT alleles
systematically closer to read edges than REF alleles?"

**Statistical method**: Same Z/√N effect-size normalization as MQCD, applied
to folded read positions instead of mapping qualities.

**Computation**:

1. During minimap2 re-alignment in the genotyper, walk the best-allele CIGAR
   to find the query position corresponding to the variant's haplotype position.
2. Compute folded position: `min(qpos/read_length, 1 - qpos/read_length)`.
   Result: 0.0 = read edge, 0.5 = read center.
3. Group folded positions by allele (REF vs ALT).
4. Apply Mann-Whitney U test → Z/√N effect size, computed per-sample.

**Value range**: [−2, +2] typically, or `0.0` (untestable)

**Coverage stability**: The effect size is coverage-invariant. A consistent
edge bias produces the same RPCD at any depth from 20× to 2000×.

**Interpretation**:

| RPCD Range | Meaning |
|:-----------|:--------|
| ≈ 0.0 | Uniform read position distribution — expected for true variants. Also returned when untestable. |
| −0.2 to −0.5 | Moderate edge bias — ALT allele somewhat closer to read ends |
| < −0.5 | Strong edge bias — likely alignment artifact or 3' error cascade |

**Manual interpretation tip**: RPCD should be interpreted jointly with BQCD.
Many artifacts produce both read-edge clustering (negative RPCD) and low
ALT base quality (negative BQCD) because base quality degrades near read
ends. If RPCD < −0.3 and BQCD < −0.3, the variant is very likely an artifact.
True variants may show mild negative RPCD in isolation (e.g., near an indel
that shifts alignment) without the accompanying BQCD depression.

---

## Base Quality Cohen's D (`BQCD`)

**Purpose**: Detect chemistry-driven sequencing artifacts where the ALT allele
is systematically called with lower base confidence than the REF allele.

**Key artifact**: 8-oxoguanine (8-oxoG) oxidation is the dominant source of
G→T / C→A errors in Illumina sequencing. Oxidized guanine mispairs with
adenine, producing consistent low-quality G→T miscalls. The miscalled base
has characteristically low Phred scores that are detectable by this test.

**Statistical method**: Z/√N effect-size normalization on per-read base
qualities at the variant position (REF vs ALT groups).

**Computation**:

1. During genotyping, the representative base quality at the variant position
   is recorded per read (minimum across the variant region for indels).
2. Base qualities are grouped by allele (REF vs ALT), combining forward and
   reverse strand observations.
3. Apply Mann-Whitney U test → Z/√N effect size, computed per-sample.

**Value range**: [−2, +2] typically, or `0.0` (untestable)

**Coverage stability**: Coverage-invariant. A consistent 5-unit ALT quality
depression produces the same BQCD at any depth.

**Interpretation**:

| BQCD Range | Meaning |
|:-----------|:--------|
| ≈ 0.0 | No systematic quality difference — expected for true variants. Also returned when untestable. |
| −0.2 to −0.5 | Moderate ALT quality depression — inspect for oxidation or deamination |
| < −0.5 | Strong signal — likely chemistry artifact (8-oxoG, FFPE deamination) |

**Manual interpretation tip**: For targeted 8-oxoG detection, examine BQCD
jointly with the variant allele: G→T and C→A substitutions with BQCD < −0.3
are the classic oxidation signature. Other mutation types with negative BQCD
may indicate FFPE deamination (C→T/U) or other library damage.

---

## Allele-Specific Mismatch Delta (`ASMD`)

**Purpose**: Detect chimeric reads and paralogous mismapping. True variants
should be the only difference between a read and the reference. If ALT-supporting
reads carry excess random mismatches while REF reads are clean, the ALT allele
is likely a misaligned chimera or paralog.

**Computation**:

1. During minimap2 re-alignment, every read is aligned to the **REF haplotype**
   (haplotype index 0), regardless of which allele it is ultimately assigned to.
2. The edit distance (NM) is computed from the REF alignment CIGAR per SAM spec:
   mismatches (under M ops) + insertion bases + deletion bases. Soft clips,
   hard clips, and reference skips are excluded.
3. NM values are grouped by allele assignment (REF vs ALT).
4. `ASMD = mean(ALT NM) − mean(REF NM)`, computed per-sample.

**Why NM cancels the variant's own contribution**: Both REF and ALT reads are
aligned to the same REF haplotype. The variant itself is an edit against the
reference for all reads equally. ASMD therefore isolates *excess noise* beyond
the shared baseline.

**Value range**: (−∞, +∞), typically [0, 20]

**Coverage stability**: Mean edit distance converges quickly. ASMD is stable
above 10× per allele. At extreme coverages (1000×+), the mean becomes very
precise, but the expected value remains the same.

**Interpretation**:

| ASMD Range | Meaning |
|:-----------|:--------|
| ≈ 0 | ALT and REF reads have similar edit distance — true variant expected |
| 1–3 | Mild excess noise — possibly repetitive region or low-quality library |
| > 5 | Strong signal — ALT reads carry many extra mismatches, likely misaligned from a paralogous locus or chimeric junction |
| 0.0 (empty) | Either REF or ALT group had no reads |

**Manual interpretation tip**: ASMD should be interpreted together with MQCD.
Paralogous mismapping usually shows both elevated ASMD (excess mismatches in
mismapped ALT reads) and depressed MQCD (low ALT mapping quality). True
variants in repetitive regions may show mild ASMD elevation (1–2) without
MQCD depression.

---

## Site Depth Fold Change (`SDFC`)

**Purpose**: Detect collapsed paralogous mappings and abnormal depth at the
variant site relative to the local background.

**Computation**:

1. During `ProcessWindow`, the window mean coverage is computed from all
   sampled BAM reads across all samples:
   `WindowCov = Σ(sampled_bases) / window_length`.
2. For each variant, `SDFC = DP / WindowCov`, where DP is the total read
   depth at the variant site (sum across alleles).

**Why window mean coverage**: The window (≥ 1000 bp, enforced minimum) averages
read depth across hundreds of positions, providing a stable local background
estimate that is immune to variant density non-uniformity and single-position
outliers. This is fundamentally more robust than variant-DP-based approaches
(e.g., EMA of nearby variant depths) which are sparse and outlier-sensitive.

**Value range**: [0, ∞)

**Coverage stability**: Inherently coverage-invariant — SDFC is a ratio (site
depth / window depth). A 2× depth spike produces SDFC ≈ 2.0 at any overall
coverage level.

**Interpretation**:

| SDFC Range | Meaning |
|:-----------|:--------|
| 0.8–1.2 | Normal depth — variant site matches local background |
| > 2.0 | Elevated depth — possible collapsed paralog or segmental duplication mapped to one locus |
| > 5.0 | Extreme — strong paralogous collapse signal |
| < 0.5 | Depleted depth — possible allelic dropout, mapping hole, or deletion |

**Manual interpretation tip**: Use SDFC together with MQCD for a comprehensive
paralog detection strategy: collapsed paralogs produce both elevated SDFC (excess
reads mapped to one site) and depressed MQCD (the paralogous reads have lower
MAPQ). A variant at a site with SDFC > 2.0 and MQCD < −0.3 is very likely a
false positive.

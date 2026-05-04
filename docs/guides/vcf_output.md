# VCF Output Reference

Lancet2 produces bgzipped VCF output following the VCF 4.5 specification.
This page is the authoritative reference for every INFO and FORMAT field
emitted, including expected value ranges and interpretation guidelines.
Detailed computation methods are linked in separate guides.

---

## Operating Modes

Lancet2 operates in two modes that affect which INFO fields are present:

| Mode | Inputs | INFO State Tags | QUAL Source |
|:-----|:-------|:----------------|:------------|
| **Somatic** (case-control) | `--normal` + (`--tumor` and/or `--sample` with role `case`) | `SHARED`, `CTRL`, `CASE` | Somatic log odds ratio (SOLOR) |
| **Germline** (normal-only) | `--normal` BAM/CRAM only | *(none — state is UNKNOWN)* | Ref-hom PL (Dirichlet-Multinomial model, naturally asymptotes) |

In somatic mode, state tags classify each variant by ALT allele
presence across sample roles. In germline mode, state classification
is not possible, so no state flags appear in the INFO field.

---

## INFO Fields

### Core Fields (always present)

| Field | VCF Type | Description | Range | Interpretation |
|:------|:---------|:------------|:------|:---------------|
| `MULTIALLELIC` | Flag | Multiallelic indicator | — | Present if the record contains more than one ALT allele at this locus. |
| `TYPE` | String | Variant type | `SNV`, `INS`, `DEL`, `MNP` | The class of variant event detected during local assembly. Comma-separated array for multiallelic variants. |
| `LENGTH` | Integer | Variant length in base pairs | ≥ 1 | For SNVs: always 1. For INDELs: the number of inserted or deleted bases. Comma-separated array for multiallelic variants. |

### Somatic State Flags (case-control mode only)

| Field | VCF Type | Description | Interpretation |
|:------|:---------|:------------|:---------------|
| `SHARED` | Flag | ALT allele seen in both case and control | Likely germline variant or LOH event |
| `CTRL` | Flag | ALT allele seen only in control (normal) sample(s) | Possible loss-of-heterozygosity in case |
| `CASE` | Flag | ALT allele seen only in case (tumor) sample(s) | Candidate somatic variant |

### Complexity Annotations
Graph topology and sequence composition metrics computed during assembly.
Coverage-invariant by design — suitable as direct features for variant
quality models without depth normalization.

| Field | VCF Type | Description |
|:------|:---------|:------------|
| `GRAPH_CX` | String (3 values) | Graph structural complexity (GEI, TipToPathCovRatio, MaxSingleDirDegree). Topology-derived; coverage-stable above 20×. See [Graph Complexity](graph_complexity.md). |
| `SEQ_CX` | String (11 values) | Sequence complexity (ContextHRun, ContextEntropy, ContextFlankLQ, ContextHaplotypeLQ, DeltaHRun, DeltaEntropy, DeltaFlankLQ, TrAffinity, TrPurity, TrPeriod, IsStutterIndel). Perfectly coverage-invariant (sequence-only). See [Sequence Complexity](sequence_complexity.md). |

---

## Sample Column Order

Sample FORMAT columns appear in **deterministic sorted order**: first by role (control before case), then by SM read group tag name (lexicographic). This order is independent of the order in which `--normal`, `--tumor`, or `--sample` flags appear on the command line. Multiple BAM/CRAM files sharing the same SM tag and role are treated as one logical sample.

---

## FORMAT Fields (per-sample)

### Coverage Stability Overview

Lancet2's FORMAT fields are designed with ML generalization in mind.
Most derived metrics are **coverage-invariant**: they produce the same
value for the same biological signal regardless of sequencing depth
(20×–2000×). Raw count fields (AD, ADF, ADR, DP) intentionally scale
with coverage — they are informational and should not be used as direct
ML features. The normalized versions (NPBQ, PRAD, PANG, SDFC, SB, RPCD,
BQCD, MQCD) are the coverage-stable alternatives.

| Stability | Fields |
|:----------|:-------|
| **Perfectly invariant** | PANG, SCA, SB, RPCD, BQCD, MQCD, FSSE, HSE |
| **Bounded, near-invariant** | NPBQ (~30), ASMD, SDFC, PRAD ([0, 3.5]), AHDD, PDCV |
| **Bounded by ceiling** | RMQ ([0, 60]), GQ ([0, 99]) |
| **Raw counts (not for ML)** | GT, AD, ADF, ADR, DP, PL |

The FORMAT column header is: `GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ`

| Field | Number | Type | Description | Range | Interpretation |
|:------|:-------|:-----|:------------|:------|:---------------|
| `GT` | 1 | String | Genotype | VCF notation | Most likely diploid genotype (e.g., `0/1`, `1/2`). |
| `AD` | R | Integer | Allele depth | [0, ∞) | Deduplicated reads per allele (REF, ALT1, ALT2, …). |
| `ADF` | R | Integer | Forward strand depth | [0, ∞) | Forward-strand reads per allele. |
| `ADR` | R | Integer | Reverse strand depth | [0, ∞) | Reverse-strand reads per allele. |
| `DP` | 1 | Integer | Total depth | [0, ∞) | Total read coverage at this site in this sample. **Raw count — scales linearly with coverage. Use SDFC or PRAD instead for ML.** |
| `RMQ` | R | Float | RMS mapping quality | [0, 60] | Root-mean-square MAPQ per allele. Low = ambiguous alignment. Bounded by MAPQ ceiling (60); effectively coverage-stable. |
| `NPBQ` | R | Float | Normalized posterior base quality | [0, ~40] | Per-read quality contribution (raw PBQ / allele depth). Coverage-invariant: ~30 for Q30 reads at any depth. |
| `SB` | 1 | Float | Strand bias log odds ratio | [−4, +4] | 0 = balanced. \|SB\| > 1 = moderate bias. Sign: +/− = ALT enriched on fwd/rev strand. |
| `SCA` | 1 | Float | Soft Clip Asymmetry | [−1.0, 1.0] | ALT minus REF soft-clip fraction. See [details](alignment_annotations.md#soft-clip-asymmetry-sca). |
| `FLD` | 1 | Float | Fragment Length Delta | (−∞, +∞) | mean ALT isize − mean REF isize. Signed: negative = ALT fragments shorter. `.` if either group has no proper pairs (untestable). See [details](alignment_annotations.md#fragment-length-delta-fld). |
| `RPCD` | 1 | Float | Read Position Cohen's D | [−2, +2] | Effect size on folded read positions (ALT vs REF). Negative = ALT at read edges. `.` if either group is empty (untestable). See [details](alignment_annotations.md#read-position-cohens-d-rpcd). |
| `BQCD` | 1 | Float | Base Quality Cohen's D | [−2, +2] | Effect size on base qualities (ALT vs REF). Negative = ALT low quality. `.` if either group is empty (untestable). See [details](alignment_annotations.md#base-quality-cohens-d-bqcd). |
| `MQCD` | 1 | Float | Mapping Quality Cohen's D | [−2, +2] | Effect size on MAPQ (ALT vs REF). Negative = ALT mismapped. `.` if either group is empty (untestable). See [details](alignment_annotations.md#mapping-quality-cohens-d-mqcd). |
| `ASMD` | 1 | Float | Allele Mismatch Delta | (−∞, +∞) | mean(ALT NM) − mean(REF NM) − variant_length. Subtracts the variant's own edit distance so only excess noise remains. `.` if either group is empty (untestable). See [details](alignment_annotations.md#allele-specific-mismatch-delta-asmd). |
| `SDFC` | 1 | Float | Site Depth Fold Change | [0, ∞) | Sample DP / per-sample window mean coverage. >2 = possible paralog collapse. See [details](alignment_annotations.md#site-depth-fold-change-sdfc). |
| `PRAD` | 1 | Float | Polar Radius | [0, ~3.5] | log10(1 + sqrt(AD_Ref² + AD_Alt²)). Coverage-invariant signal magnitude. See [details](polar_features.md#polar-radius-prad). |
| `PANG` | 1 | Float | Polar Angle | [0, π/2] | atan2(AD_Alt, AD_Ref) in radians. 0 = hom REF, π/4 = het, π/2 = hom ALT. See [details](polar_features.md#polar-angle-pang). |
| `CMLOD` | A | Float | Continuous Mixture LOD | [0, ∞) | Per-ALT base-quality-weighted log-odds score. Integrates exact per-read base qualities into a K-dimensional mixture model comparing the MLE frequency vs. a null hypothesis. See [details](variant_discovery_genotyping.md#continuous-mixture-log-odds-cmlod). |
| `FSSE` | 1 | Float | Fragment Start Shannon Entropy | [0, 1] | Spatial diversity of ALT read start positions, normalized by log₂(min(N, 20)). Detects PCR duplicates that survive MarkDuplicates. `.` if < 3 ALT reads. See [details](alignment_annotations.md#fragment-start-shannon-entropy-fsse). |
| `AHDD` | 1 | Float | ALT-Haplotype Discordance Delta | (−∞, +∞) | mean(ALT NM vs own haplotype) − mean(REF NM vs REF). Detects assembly hallucinations where the ALT path does not actually match its supporting reads. `.` if either group is empty. See [details](alignment_annotations.md#alt-haplotype-discordance-delta-ahdd). |
| `HSE` | 1 | Float | Haplotype Segregation Entropy | [0, 1] | How concentrated ALT reads are on a single SPOA path, normalized by log₂(total_haplotypes). Near 0 = concentrated (true variant). Near 1 = scattered (noise). `.` if < 3 ALT reads or single haplotype. See [details](alignment_annotations.md#haplotype-segregation-entropy-hse). |
| `PDCV` | 1 | Float | Path Depth Coefficient of Variation | [0, ∞) | K-mer coverage uniformity along the ALT de Bruijn graph path. Max across ALT alleles. High values signal chimeric junctions with uneven support. `.` if path has < 2 nodes. See [details](alignment_annotations.md#path-depth-coefficient-of-variation-pdcv). |
| `PL` | G | Integer | Phred-scaled genotype likelihoods | [0, ∞) | Lower = more likely genotype. Computed via the Dirichlet-Multinomial model — PLs plateau at high depth due to overdispersion. PL[0]=RR, PL[1]=RA, PL[2]=AA. |
| `GQ` | 1 | Integer | Genotype quality | [0, 99] | Second-lowest Dirichlet-Multinomial PL, capped at 99. Reaches cap quickly at ≥30×; effectively coverage-stable above this threshold. |

### QUAL Column

The VCF QUAL field depends on operating mode:

- **Case-Control mode**: Somatic log odds ratio (SOLOR) comparing ALT
  enrichment between case and control:
  `QUAL = ln(((case_alt+1)(ctrl_ref+1)) / ((case_ref+1)(ctrl_alt+1)))`. Uses
  Haldane correction (+1) for zero-count robustness. Coverage-invariant:
  measures effect size, not statistical significance. Range: [0, ~10].
  `QUAL > 4` = strong somatic evidence at any coverage.

- **Control-Only mode**: Ref-hom PL directly (PL[0/0] from the Dirichlet-Multinomial model).
  The Dirichlet-Multinomial overdispersion parameter causes PLs to asymptote
  naturally at high depth — no artificial capping or depth-normalization is
  needed. Takes the maximum across samples (strongest evidence wins).

---

## Example VCF Records

### Case-Control Mode

```
chr1  12345  .  A  AT  4.85  .  CASE;TYPE=INS;LENGTH=1;GRAPH_CX=0.85,0.01,3;SEQ_CX=20,0.50,0.69,0.22,3,−0.10,0.05,1.00,0.95,1,1  GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ  0/1:20,15:10,8:10,7:35:55.0,52.3:30.2,28.5:0.150:0.1200:15.3:-0.3200:-0.1500:-0.2800:0.450:1.00:1.3979:0.6435:12.5432:0.8521:0.350:0.1200:0.4500:0,42,180:42  0/0:30,0:15,0:15,0:30:58.1,0.0:29.8,0.0:0.000:0.0000:.:.:.:.:.:0.86:1.4914:0.0000:0.0000:.:.:.:.:0,0,270:99
```

**Per-sample annotation breakdown**:

| Field | Case (sample 1) | Control (sample 2) | Interpretation |
|:------|:-----------------|:-------------------|:---------------|
| `NPBQ` | 30.2, 28.5 | 29.8, 0.0 | Good per-read quality (~Q30) for both alleles in case. |
| `SB` | 0.150 | 0.000 | Minimal strand bias in case. Control: no ALT reads, Haldane correction gives 0.0. |
| `MQCD` | −0.2800 | `.` | Slight ALT MAPQ depression in case (within normal range). Control: untestable (no ALT reads → `.`). |
| `PRAD` | 1.3979 | 1.4914 | log10(1+25) ≈ 1.40, log10(1+30) ≈ 1.49. Both in moderate-confidence range (~1.2–1.8). |
| `PANG` | 0.6435 | 0.0000 | Case ≈37° (37% VAF). Control 0° — no ALT, consistent with somatic variant. |

### Control-Only Mode

```
chr1  12345  .  A  AT  3.00  .  TYPE=INS;LENGTH=1;GRAPH_CX=0.42,0.00,2;SEQ_CX=8,1.20,0.31,0.18,1,0.05,0.02,0.80,0.90,1,0  GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:CMLOD:FSSE:AHDD:HSE:PDCV:PL:GQ  0/1:25,10:13,5:12,5:35:56.2,50.8:29.5,27.8:0.095:0.0500:8.2:-0.1500:-0.0800:-0.2100:0.320:1.00:1.4314:0.3805:8.2100:0.7200:0.120:0.0800:0.3200:0,30,150:30
```

### Native Multiallelic Mode

```
chr1  12345  .  A  AT,G  4.85  .  SHARED;MULTIALLELIC;TYPE=INS,SNV;LENGTH=1,1;GRAPH_CX=1.20,0.03,4;SEQ_CX=15,0.85,0.45,0.20,2,−0.05,0.03,0.95,0.88,2,0  GT:AD:ADF:ADR:DP:RMQ:NPBQ  1/2:10,15,5:5,8,3:5,7,2:30:55.0,52.3,60.0:30.2,28.5,35.0  0/1:25,5,0:15,3,0:10,2,0:30:58.1,50.0,0.0:29.8,25.0,0.0
```

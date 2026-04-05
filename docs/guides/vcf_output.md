# VCF Output Reference

Lancet2 produces bgzipped VCF output following the VCF 4.3 specification.
This page is the authoritative reference for every INFO and FORMAT field
emitted, including expected value ranges and interpretation guidelines.
Detailed computation methods are linked in separate guides.

---

## Operating Modes

Lancet2 operates in two modes that affect which INFO fields are present:

| Mode | Inputs | INFO State Tags | QUAL Source |
|:-----|:-------|:----------------|:------------|
| **Tumor-Normal** (somatic) | ≥ 1 tumor AND ≥ 1 normal BAM/CRAM | `SHARED`, `NORMAL`, `TUMOR` | Somatic Fisher score (right-tailed) |
| **Normal-Only** (germline) | Normal BAM/CRAM only | *(none — state is UNKNOWN)* | REF-homozygous PL |

In tumor-normal mode, state tags classify each variant by ALT allele
presence across sample types. In normal-only mode, state classification
is not possible, so no state flags appear in the INFO field.

---

## INFO Fields

### Core Fields (always present)

| Field | VCF Type | Description | Range | Interpretation |
|:------|:---------|:------------|:------|:---------------|
| `TYPE` | String | Variant type | `SNV`, `INS`, `DEL`, `MNP` | The class of variant event detected during local assembly. |
| `LENGTH` | Integer | Variant length in base pairs | ≥ 1 | For SNVs: always 1. For INDELs: the number of inserted or deleted bases. |

### Somatic State Flags (tumor-normal mode only)

| Field | VCF Type | Description | Interpretation |
|:------|:---------|:------------|:---------------|
| `SHARED` | Flag | ALT allele seen in both tumor and normal | Likely germline variant or LOH event |
| `NORMAL` | Flag | ALT allele seen only in normal sample(s) | Possible loss-of-heterozygosity in tumor |
| `TUMOR` | Flag | ALT allele seen only in tumor sample(s) | Candidate somatic variant |

### Optional Complexity Annotations

These fields require explicit CLI flags to enable. They are designed for
downstream machine learning–based variant filtering.

| Field | VCF Type | CLI Flag | Description |
|:------|:---------|:---------|:------------|
| `GRAPH_CX` | String (3 values) | `--enable-graph-complexity-features` | Graph structural complexity (GEI, TipToPathCovRatio, MaxSingleDirDegree). See [Graph Complexity](graph_complexity.md). |
| `SEQ_CX` | String (11 values) | `--enable-sequence-complexity-features` | Sequence complexity (ContextHRun, ContextEntropy, ContextFlankLQ, ContextHaplotypeLQ, DeltaHRun, DeltaEntropy, DeltaFlankLQ, TrAffinity, TrPurity, TrPeriod, IsStutterIndel). See [Sequence Complexity](sequence_complexity.md). |

---

## FORMAT Fields (per-sample)

The FORMAT column header is: `GT:AD:ADF:ADR:DP:RMQ:PBQ:SB:SCA:FLD:RPRS:BQRS:MQRS:ASMD:SDFC:PRAD:PANG:PL:GQ`

| Field | Number | Type | Description | Range | Interpretation |
|:------|:-------|:-----|:------------|:------|:---------------|
| `GT` | 1 | String | Genotype | VCF notation | Most likely diploid genotype (e.g., `0/1`, `1/2`). |
| `AD` | R | Integer | Allele depth | [0, ∞) | Deduplicated reads per allele (REF, ALT1, ALT2, …). |
| `ADF` | R | Integer | Forward strand depth | [0, ∞) | Forward-strand reads per allele. |
| `ADR` | R | Integer | Reverse strand depth | [0, ∞) | Reverse-strand reads per allele. |
| `DP` | 1 | Integer | Total depth | [0, ∞) | Total read coverage at this site in this sample. |
| `RMQ` | R | Float | RMS mapping quality | [0, 60] | Root-mean-square MAPQ per allele. Low = ambiguous alignment. |
| `PBQ` | R | Float | Posterior base quality | [0, ~60] | Bayesian per-allele base confidence. Higher = more confident. |
| `SB` | 1 | Float | Phred-scaled strand bias | [0, ∞) | 0 = no bias. > 10 = moderate. > 20 = strong. |
| `SCA` | 1 | Float | Soft Clip Asymmetry | [−1.0, 1.0] | ALT minus REF soft-clip fraction. Positive = ALT disproportionately clipped. See [details](alignment_annotations.md#soft-clip-asymmetry-sca). |
| `FLD` | 1 | Float | Fragment Length Delta | [0, ∞) | \|mean ALT isize − mean REF isize\|. Large values flag chimeric artifacts. See [details](alignment_annotations.md#fragment-length-delta-fld). |
| `RPRS` | 1 | Float | Read Position Rank Sum | (−∞, +∞) or 100.0 | Mann-Whitney U Z-score on folded read positions (ALT vs REF). Negative = ALT at read edges. See [details](alignment_annotations.md#read-position-rank-sum-z-score-rprs). |
| `BQRS` | 1 | Float | Base Quality Rank Sum | (−∞, +∞) or 100.0 | Mann-Whitney U Z-score on base qualities (ALT vs REF). Negative = ALT low quality. See [details](alignment_annotations.md#base-quality-rank-sum-z-score-bqrs). |
| `MQRS` | 1 | Float | Mapping Quality Rank Sum | (−∞, +∞) or 100.0 | Mann-Whitney U Z-score (ALT vs REF MAPQ). Negative = ALT mismapped. 100.0 = test not applicable. See [details](alignment_annotations.md#mapping-quality-rank-sum-z-score-mqrs). |
| `ASMD` | 1 | Float | Allele Mismatch Delta | (−∞, +∞) | mean(ALT NM) − mean(REF NM) from REF haplotype alignment. Positive = ALT excess noise. See [details](alignment_annotations.md#allele-specific-mismatch-delta-asmd). |
| `SDFC` | 1 | Float | Site Depth Fold Change | [0, ∞) | DP / window mean coverage. >2 = possible paralog collapse. <0.5 = mapping hole. See [details](alignment_annotations.md#site-depth-fold-change-sdfc). |
| `PRAD` | 1 | Float | Polar Radius | [0, ∞) | sqrt(AD_Ref² + AD_Alt²). Signal magnitude orthogonal to allele fraction. See [details](polar_features.md#polar-radius-prad). |
| `PANG` | 1 | Float | Polar Angle | [0, π/2] | atan2(AD_Alt, AD_Ref) in radians. 0 = hom REF, π/4 = het, π/2 = hom ALT. See [details](polar_features.md#polar-angle-pang). |
| `PL` | G | Integer | Phred-scaled genotype likelihoods | [0, ∞) | Lower = more likely genotype. PL[0]=RR, PL[1]=RA, PL[2]=AA. |
| `GQ` | 1 | Integer | Genotype quality | [0, 99] | Second-lowest PL, capped at 99. 0 = tied PLs. |

### QUAL Column

The VCF QUAL field depends on operating mode:

- **Tumor-Normal mode**: Somatic Fisher score — Phred-scaled right-tailed
  Fisher's exact test p-value from a 2×2 table of (tumor ALT, tumor REF)
  vs (mean normal ALT, mean normal REF). Uses the right-tail (one-sided)
  because somatic calling tests the directional hypothesis: "tumor has
  higher ALT enrichment than normal."

- **Normal-Only mode**: The REF-homozygous PL value (`PL[0]`), representing
  confidence that the site differs from reference.

---

## Example VCF Records

### Tumor-Normal Mode

```
chr1  12345  .  A  AT  42.50  .  TUMOR;TYPE=INS;LENGTH=1  GT:AD:ADF:ADR:DP:RMQ:PBQ:SB:SCA:FLD:MQRS:PRAD:PANG:PL:GQ  0/1:20,15:10,8:10,7:35:55.0,52.3:38.5,36.1:4.2:0.1200:15.3:-0.450:25.0:0.6435:0,42,180:42  0/0:30,0:15,0:15,0:30:58.1,0.0:40.2,0.0:0.0:0.0000:0.0:100.000:30.0:0.0000:0,0,270:99
```

**Per-sample annotation breakdown**:

| Field | Tumor (sample 1) | Normal (sample 2) | Interpretation |
|:------|:-----------------|:-------------------|:---------------|
| `SCA` | 0.1200 | 0.0000 | Tumor ALT reads disproportionately soft-clipped; normal has no ALT reads. |
| `FLD` | 15.3 | 0.0 | Moderate insert size discrepancy in tumor; normal has no ALT evidence. |
| `MQRS` | −0.450 | 100.000 | Slight ALT MAPQ depression in tumor (within normal range). Normal: sentinel (no ALT reads). |
| `PRAD` | 25.0 | 30.0 | Both samples have moderate depth. Normal is slightly deeper. |
| `PANG` | 0.6435 | 0.0000 | Tumor ≈37° (37% VAF). Normal 0° — no ALT, consistent with somatic variant. |

### Normal-Only Mode

```
chr1  12345  .  A  AT  30.00  .  TYPE=INS;LENGTH=1  GT:AD:ADF:ADR:DP:RMQ:PBQ:SB:SCA:FLD:MQRS:PRAD:PANG:PL:GQ  0/1:25,10:13,5:12,5:35:56.2,50.8:39.0,35.5:2.1:0.0500:8.2:-0.230:26.9:0.3805:0,30,150:30
```

### With All Annotations Enabled

```
chr1  12345  .  A  AT  42.50  .  TUMOR;TYPE=INS;LENGTH=1;GRAPH_CX=0.85,0.01,3;SEQ_CX=20,0.50,0.69,0.22,3,−0.10,0.05,1.00,0.95,1,1  GT:AD:ADF:ADR:DP:RMQ:PBQ:SB:SCA:FLD:MQRS:PRAD:PANG:PL:GQ  0/1:20,15:10,8:10,7:35:55.0,52.3:38.5,36.1:4.2:0.1200:15.3:-0.450:25.0:0.6435:0,42,180:42
```

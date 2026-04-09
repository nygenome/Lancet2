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
| **Tumor-Normal** (somatic) | ≥ 1 tumor AND ≥ 1 normal BAM/CRAM | `SHARED`, `NORMAL`, `TUMOR` | Somatic log odds ratio (SOLOR) |
| **Normal-Only** (germline) | Normal BAM/CRAM only | *(none — state is UNKNOWN)* | Per-read evidence: min(PL[0]/DP, 10) |

In tumor-normal mode, state tags classify each variant by ALT allele
presence across sample types. In normal-only mode, state classification
is not possible, so no state flags appear in the INFO field.

---

## INFO Fields

### Core Fields (always present)

| Field | VCF Type | Description | Range | Interpretation |
|:------|:---------|:------------|:------|:---------------|
| `MULTIALLELIC` | Flag | Multiallelic indicator | — | Present if the record natively aggregates more than one ALT structural allele evaluated at the exact same locus constraints. |
| `TYPE` | String | Variant type | `SNV`, `INS`, `DEL`, `MNP` | The class of variant event detected during local assembly. Comma-separated array for multiallelic variants. |
| `LENGTH` | Integer | Variant length in base pairs | ≥ 1 | For SNVs: always 1. For INDELs: the number of inserted or deleted bases. Comma-separated array for multiallelic variants. |

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
| `GRAPH_CX` | String (3 values) | `--enable-graph-complexity-features` | Graph structural complexity (GEI, TipToPathCovRatio, MaxSingleDirDegree). Topology-derived; coverage-stable above 20×. See [Graph Complexity](graph_complexity.md). |
| `SEQ_CX` | String (11 values) | `--enable-sequence-complexity-features` | Sequence complexity (ContextHRun, ContextEntropy, ContextFlankLQ, ContextHaplotypeLQ, DeltaHRun, DeltaEntropy, DeltaFlankLQ, TrAffinity, TrPurity, TrPeriod, IsStutterIndel). Perfectly coverage-invariant (sequence-only). See [Sequence Complexity](sequence_complexity.md). |

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
| **Perfectly invariant** | PANG, SCA, SB, RPCD, BQCD, MQCD |
| **Bounded, near-invariant** | NPBQ (~30), ASMD, SDFC, PRAD ([0, 3.5]) |
| **Bounded by ceiling** | RMQ ([0, 60]), GQ ([0, 99]) |
| **Raw counts (not for ML)** | GT, AD, ADF, ADR, DP, PL |

The FORMAT column header is: `GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ`

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
| `FLD` | 1 | Float | Fragment Length Delta | [0, ∞) | \|mean ALT isize − mean REF isize\|. See [details](alignment_annotations.md#fragment-length-delta-fld). |
| `RPCD` | 1 | Float | Read Position Cohen's D | [−2, +2] | Effect size on folded read positions (ALT vs REF). Negative = ALT at read edges. 0 = untestable. See [details](alignment_annotations.md#read-position-cohens-d-rpcd). |
| `BQCD` | 1 | Float | Base Quality Cohen's D | [−2, +2] | Effect size on base qualities (ALT vs REF). Negative = ALT low quality. 0 = untestable. See [details](alignment_annotations.md#base-quality-cohens-d-bqcd). |
| `MQCD` | 1 | Float | Mapping Quality Cohen's D | [−2, +2] | Effect size on MAPQ (ALT vs REF). Negative = ALT mismapped. 0 = untestable. See [details](alignment_annotations.md#mapping-quality-cohens-d-mqcd). |
| `ASMD` | 1 | Float | Allele Mismatch Delta | (−∞, +∞) | mean(ALT NM) − mean(REF NM). Positive = ALT excess noise. See [details](alignment_annotations.md#allele-specific-mismatch-delta-asmd). |
| `SDFC` | 1 | Float | Site Depth Fold Change | [0, ∞) | DP / window mean coverage. >2 = possible paralog collapse. See [details](alignment_annotations.md#site-depth-fold-change-sdfc). |
| `PRAD` | 1 | Float | Polar Radius | [0, ~3.5] | log10(1 + sqrt(AD_Ref² + AD_Alt²)). Coverage-invariant signal magnitude. See [details](polar_features.md#polar-radius-prad). |
| `PANG` | 1 | Float | Polar Angle | [0, π/2] | atan2(AD_Alt, AD_Ref) in radians. 0 = hom REF, π/4 = het, π/2 = hom ALT. See [details](polar_features.md#polar-angle-pang). |
| `PL` | G | Integer | Phred-scaled genotype likelihoods | [0, ∞) | Lower = more likely genotype. PL[0]=RR, PL[1]=RA, PL[2]=AA. **Scales with coverage** — the QUAL column provides coverage-normalized versions. |
| `GQ` | 1 | Integer | Genotype quality | [0, 99] | Second-lowest PL, capped at 99. Reaches cap quickly at ≥30×; effectively coverage-stable above this threshold. |

### QUAL Column

The VCF QUAL field depends on operating mode:

- **Tumor-Normal mode**: Somatic log odds ratio (SOLOR) comparing ALT
  enrichment between tumor and normal:
  `QUAL = ln(((tmr_alt+1)(nml_ref+1)) / ((tmr_ref+1)(nml_alt+1)))`. Uses
  Haldane correction (+1) for zero-count robustness. Coverage-invariant:
  measures effect size, not statistical significance. Range: [0, ~10].
  `QUAL > 4` = strong somatic evidence at any coverage.

- **Normal-Only mode**: Per-read evidence against hom-REF, computed as
  `min(max(PL[0]/DP per sample), 10.0)`. Coverage-invariant: a clean het
  produces QUAL ≈ 3.0 at any depth. Takes the maximum across samples
  (strongest evidence wins).

---

## Example VCF Records

### Tumor-Normal Mode

```
chr1  12345  .  A  AT  4.85  .  TUMOR;TYPE=INS;LENGTH=1  GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ  0/1:20,15:10,8:10,7:35:55.0,52.3:30.2,28.5:0.150:0.1200:15.3:-0.3200:-0.1500:-0.2800:0.450:1.00:1.3979:0.6435:0,42,180:42  0/0:30,0:15,0:15,0:30:58.1,0.0:29.8,0.0:0.000:0.0000:0.0:0.0000:0.0000:0.0000:0.000:0.86:1.4914:0.0000:0,0,270:99
```

**Per-sample annotation breakdown**:

| Field | Tumor (sample 1) | Normal (sample 2) | Interpretation |
|:------|:-----------------|:-------------------|:---------------|
| `NPBQ` | 30.2, 28.5 | 29.8, 0.0 | Good per-read quality (~Q30) for both alleles in tumor. |
| `SB` | 0.150 | 0.000 | Minimal strand bias in tumor. Normal: no ALT reads, Haldane correction gives 0.0. |
| `MQCD` | −0.2800 | 0.0000 | Slight ALT MAPQ depression in tumor (within normal range). Normal: untestable (no ALT). |
| `PRAD` | 1.3979 | 1.4914 | log10(1+25) ≈ 1.40, log10(1+30) ≈ 1.49. Both in moderate-confidence range (~1.2–1.8). |
| `PANG` | 0.6435 | 0.0000 | Tumor ≈37° (37% VAF). Normal 0° — no ALT, consistent with somatic variant. |

### Normal-Only Mode

```
chr1  12345  .  A  AT  3.00  .  TYPE=INS;LENGTH=1  GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ  0/1:25,10:13,5:12,5:35:56.2,50.8:29.5,27.8:0.095:0.0500:8.2:-0.1500:-0.0800:-0.2100:0.320:1.00:1.4314:0.3805:0,30,150:30
```

### With All Annotations Enabled

```
chr1  12345  .  A  AT  4.85  .  TUMOR;TYPE=INS;LENGTH=1;GRAPH_CX=0.85,0.01,3;SEQ_CX=20,0.50,0.69,0.22,3,−0.10,0.05,1.00,0.95,1,1  GT:AD:ADF:ADR:DP:RMQ:NPBQ:SB:SCA:FLD:RPCD:BQCD:MQCD:ASMD:SDFC:PRAD:PANG:PL:GQ  0/1:20,15:10,8:10,7:35:55.0,52.3:30.2,28.5:0.150:0.1200:15.3:-0.3200:-0.1500:-0.2800:0.450:1.00:1.3979:0.6435:0,42,180:42
```

### Native Multiallelic Mode

```
chr1  12345  .  A  AT,G  4.85  .  SHARED;MULTIALLELIC;TYPE=INS,SNV;LENGTH=1,1  GT:AD:ADF:ADR:DP:RMQ:NPBQ  1/2:10,15,5:5,8,3:5,7,2:30:55.0,52.3,60.0:30.2,28.5,35.0  0/1:25,5,0:15,3,0:10,2,0:30:58.1,50.0,0.0:29.8,25.0,0.0
```

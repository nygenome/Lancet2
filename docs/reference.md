---
hide:
  - navigation
---

# Reference

## pipeline
The pipeline subcommand runs the entire Lancet variant calling pipeline on one (or) more region(s) of interest.
The full help text for the subcommand can be generated using the following command line

```bash
Lancet2 pipeline --help
```

### Required

#### `-r`,`--reference`
> [PATH]

Path to the reference FASTA file

#### `-o`,`--out-vcfgz`
> [PATH]

Output path to the compressed VCF file

### Datasets

#### `-n`,`--normal`
> [PATH...]

Path to one (or) more normal BAM/CRAM file(s). Required.

#### `-t`,`--tumor`
> [PATH...]

Path to one (or) more tumor BAM/CRAM file(s). Optional — when omitted, Lancet2 runs in germline-only mode.

### Regions

#### `-R`,`--region`
> [REF:[:START[-END]]...]

One (or) more regions (1-based both inclusive).

#### `-b`,`--bed-file`
> [PATH]

Path to BED file with regions to process

#### `-P`,`--padding`
> [0-1000]. Default value --> 500

Padding for both sides of all input regions

#### `-p`,`--pct-overlap`
> [10-90]. Default value --> 20

Percent overlap between consecutive windows

#### `-w`,`--window-size`
> [250-2500]. Default value --> 1000

Window size for micro-assembly and variant calling tasks

### Parameters

#### `-T`,`--num-threads`
Number of additional async worker threads. Default value --> 2

#### `-k`,`--min-kmer`
Minimum kmer length to try for micro-assembly graph nodes. Default value --> 13

#### `-K`,`--max-kmer`
Maximum kmer length to try for micro-assembly graph nodes. Default value --> 127

#### `-s`,`--kmer-step`
Kmer step size length to try for micro-assembly graph nodes. Must be one of {2, 4, 6, 8, 10}. Default value --> 6

#### `--min-anchor-cov`
Minimum coverage for micro-assembly graph anchor nodes (source/sink). Default value --> 5

#### `--min-node-cov`
Minimum coverage for nodes in the micro-assembly graph. Default value --> 2

#### `--max-sample-cov`
Maximum per sample coverage before downsampling. Default value --> 1000

### Flags

#### `--verbose`
Turn on verbose logging

#### `--extract-pairs`
Extract all useful read pairs

#### `--no-active-region`
Force assemble all windows

#### `--no-contig-check`
Skip contig check with reference

### Optional

#### `--graphs-dir`
Output directory to write per window graphs in DOT and GFA format.
Must be a non-existing directory path that will be created.

#### `--annotation-features`
Comma-separated list of optional annotation feature IDs to include in VCF output.
Supported features: `ALT_LCR`, `REF_LCR`

```bash
Lancet2 pipeline --annotation-features ALT_LCR,REF_LCR \
    --normal normal.bam --tumor tumor.bam \
    --reference ref.fasta --region "chr22" \
    --out-vcfgz output.vcf.gz
```

## VCF Output

### INFO fields

| Field      | Type    | Description |
|:-----------|:--------|:------------|
| `SHARED`   | Flag    | Variant ALT seen in both tumor & normal sample(s). Somatic mode only. |
| `NORMAL`   | Flag    | Variant ALT seen only in normal sample(s). Somatic mode only. |
| `TUMOR`    | Flag    | Variant ALT seen only in tumor sample(s). Somatic mode only. |
| `TYPE`     | String  | Variant type. Possible values: `SNV`, `INS`, `DEL`, `MNP` |
| `LENGTH`   | Integer | Variant length in base pairs |
| `ALT_LCR`  | Float   | Low-complexity scores at 50bp, 100bp flanks (k=4) and full ALT haplotype (k=7). Optional, requires `--annotation-features ALT_LCR`. |
| `REF_LCR`  | Float   | Low-complexity scores at 50bp, 100bp flanks (k=4) and full REF haplotype (k=7). Optional, requires `--annotation-features REF_LCR`. |

### FORMAT fields

| Field | Number | Type    | Description |
|:------|:-------|:--------|:------------|
| `GT`  | 1      | String  | Genotype called at the variant site |
| `AD`  | R      | Integer | Read depth per allele (REF, ALT1, ALT2, ...) |
| `ADF` | R      | Integer | Forward strand read depth per allele |
| `ADR` | R      | Integer | Reverse strand read depth per allele |
| `DP`  | 1      | Integer | Total read depth at the variant site |
| `RMQ` | R      | Float   | RMS mapping quality per allele |
| `PBQ` | R      | Float   | Posterior base quality per allele (Bayesian aggregation) |
| `SB`  | R      | Float   | Strand bias ratio per allele (fwd/total) |
| `PL`  | G      | Integer | Phred-scaled genotype likelihoods |
| `GQ`  | 1      | Integer | Genotype quality (second-lowest PL, capped at 99) |

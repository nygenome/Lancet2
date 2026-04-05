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

#### `--enable-graph-complexity-features`
Emit `GRAPH_CX` INFO tag with per-variant graph complexity metrics.
See [VCF Output Reference](guides/vcf_output.md#graph-complexity-graph_cx) for details.

#### `--enable-sequence-complexity-features`
Emit `ULTRA_*_CX`, `MICRO_*_CX`, and `MACRO_REF_CX` INFO tags with multi-scale
sequence complexity metrics (HRun, entropy, TR motifs, LongdustQ).
See [VCF Output Reference](guides/vcf_output.md#multi-scale-sequence-complexity) for details.

#### `--genome-gc-bias`
> [0.0-1.0]. Default value --> 0.41

Global genome GC fraction for LongdustQ score correction.
Set to 0.5 to disable GC correction (uniform model).
See [VCF Output Reference](guides/vcf_output.md#gc-bias-correction) for details.

```bash
Lancet2 pipeline \
    --enable-graph-complexity-features \
    --enable-sequence-complexity-features \
    --genome-gc-bias 0.41 \
    --normal normal.bam --tumor tumor.bam \
    --reference ref.fasta --region "chr22" \
    --out-vcfgz output.vcf.gz
```

## VCF Output

See the [VCF Output Format](guides/vcf_output.md) guide for complete documentation
of all INFO and FORMAT fields.


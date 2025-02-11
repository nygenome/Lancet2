---
hide:
  - navigation
---

# Reference

## pipeline
The pipeline subcommand runs the entire Lancet variant calling pipeline on one (or) more region(s) of intereset.
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

Path to one (or) more normal BAM/CRAM file(s)

#### `-t`,`--tumor`
> [PATH...]

Path to one (or) more tumor BAM/CRAM file(s)

### Regions

#### `-R`,`--region`
> [SAMTOOLS_STYLE_REGION_STRING...]

One (or) more regions (1-based both inclusive).

#### `-b`,`--bed-file`
> [PATH]

Path to BED file with regions to process

#### `-P`,`--padding`
> [0-1000]. Default value --> 500

Padding for both sides of all input regions

#### `-p`,`--pct-overlap`
> [20-80]. Default value --> 20

Percent overlap between consecutive windows

#### `-w`,`--window-size`
> [500-5000]. Default value --> 500

Window size for micro-assembly and variant calling tasks

### Parameters

#### `-T`,`--num-threads`
Number of additional async worker threads

#### `-k`,`--min-kmer`
Minimum kmer length to try for micro-assembly graph nodes

#### `-K`,`--max-kmer`
Maximum kmer length to try for micro-assembly graph nodes

#### `-s`,`--kmer-step`
Kmer step size length to try for micro-assembly graph nodes

#### `--min-anchor-cov`
Minimum coverage for micro-assembly graph anchor nodes (source/sink)

#### `--min-node-cov`
Minimum coverage for nodes in the micro-assembly graph

#### `--max-sample-cov`
Maximum per sample coverage before downsampling

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

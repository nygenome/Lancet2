---
id: cli
title: Command Line Reference
---
## pipeline
This is the subcommand that will kick off the tool. This will always follow directly after the call to the Lancet executable

```bash
lancet2 pipeline --help
lancet2 pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out_prefix
```

## General Options:
### `--help`
Bring up the standard help output. This will give examples of options the user can change to customize a run.

### `--version`
Print the version information for the build

## Required Arguments:
A standard run of Lancet will provide a tumor bam file, a normal bam file, the reference fasta, and an output path to place the outputted vcf file.
```shell
./lancet pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out.vcf
```

### `-t`, `--tumor`
Provide the path to the tumor bam file. Index for this file should also be in the same directory

### `-n`, `--normal`
Provide the path to the normal bam file. Index for this file should also be in the same directory

### `-r`, `--reference`
Provide the path to the reference fasta file. Index for this file should also be in the same directory

### `-o`, `--out-prefix`
Prefix to use for output VCF (will be bgzipped and indexed)

## Optional Arguments:
These arguments allow for more fine-tuned control of the tool. If not provided, default values will be assigned
### `--graphs-dir`
This tag allows you to define the output path for dumping serialized graphs from a run. If this option is not utilized, there will be no outputted graphs.

### Regions
These options will allow you to play around with what the tool looks at.

### `--region`
Allows the user to define what region the tool should run on. Should be of format chr:start_pos-end_pos.
For example, this will indicate that the tool will run on chromosome 2 between positions 33091000 and 33092000:
```shell
... --region 2:33091000-33092000 ...
```
If no region(s) specified, the tool will default to run on the whole genome provided

### `-b`, `--bed-file`
Path to bed file that will be used to define which region the tool will run on. Bed file should be a tab delimited file in which the first column
is the chromosome, the second column is the start position, and the third column is the end position. For example, a bed file like this specifies chromosomes
1 and 2 between positions 4 and 239410 and 892 and 1029348 respectively:
```shell
1	4	239410
2	892	1029348
```
If no region(s) specified, the tool will default to run on the whole genome provided

### `-P`, `--padding`
Padding to be applied to all input regions. By default, there will be a padding of 250 bp.

### `-w`, `--window-size`
This tag allows you to define the window size used for each microassembly task. By default, the window will be 600 bp.

### `-T`, `--num-threads`
Allows the user to specify the number of cores to be used for the tool.

### `--pct-overlap`
Allows the user to define how much overlap there should be between windows. If not specified, the tool will default to 50% overlap

### Parameters
These options allow you to define certain parameters for how the tool performs its variant calling

### `-T`, `--num-threads`
This allows you to define how many threads are used by the tool. If not specified, the tool will default to running with 1 thread.

### `-k`, `--min-kmer-length`
This allows you to define the minimum length kmers should be for graph nodes. If no length specified, the min kmer length defaults to 11 bp.

### `-K`, `--max-kmer-length`
This allows you to define the maximum length kmers should be for graph nodes. If no length specified, the max kmer length defaults to 101 bp.

### `--min-trim-qual`
This allows you to define the minimum base quality for trimming 5' and 3' read bases. If this option is not used, the minimum base quality will be defaulted to 10.

### `-q`, `--min-base-qual`
This allows you to define the minimum base quality for SNV calling. If this option is not used, the minimum base quality will be defaulted to 17.

### `-Q`, `--min-mapping-qual`
This allows you to define the minimum mapping quality required to use a read. By default, this value will be set to 15.

### `--max-rpt-mismatch`
This allows you to define the maximum number of mismatches used to detect approximate repeats. By default, the maximum number of mismatches is set to 2.

### `--max-tip-length`
This will define the maximum tip length allowed in the genereated graphs. By default, the maximum tip length is set to 11.

### `--graph-traversal-limit`
Maximum allowed tip length in the graph. By default, this value is set to 11.
Check cli.cpp, why is params->minGraphTipLength

### `--graph-traversal-limit`
Max BFS/DFS graph traversal limit. By default, this value is set to 100000.

### `--max-indel-length`
Maximum limit on the indel length to detect. Default is set to 500.

### `--min-anchor-cov`
Minimum coverage for anchor nodes (source & sink). Default is 5

### `--min-node-cov`
Minimum coverage for all nodes in the graph. Default is 1

### `--min-cov-ratio`
Minimum node by window coverage for all nodes. Default is 0.01
Node to window coverage ratio?

### `--max-window-cov`
Maximum average window coverage before downsampling. Default is 1000

### `--min-as-xs-diff`
Minimum difference between AS and XS scores (BWA-mem). Default is 5

### STR Parameters
Use these flags to deal with Short Tandem Repeats

### `--max-str-unit-len`
Maximum unit length for an STR motif. Default is 4

### `--min-str-units`
Minimum number of STR units required to report. Default is 7

### `--maxSTRDist`
Maximum distance (in bp) of variant from the STR motif. Default is 1

### Filters
These options let you apply different filters to apply to the variant caller

### `-c`, `--max-nml-alt-cnt`
Maximum ALT allele count in normal sample. Default is set at 0

### `-C`, `--min-tmr-alt-cnt`
Minimum ALT allele count in tumor sample. Default is 3

### `-v`, `--max-nml-vaf`
Maximum variant allele frequency in normal sample. Default is 0

### `-V`, `--min-tmr-vaf`
Minimum variant allele frequency in tumor sample. Default is 0.01

### `--min-nml-cov`
Minimum variant coverage in the normal sample. Default is 10

### `--min-tmr-cov`
Minimum variant coverage in the tumor sample. Default is 4

### `--max-nml-cov`
Maximum variant coverage in the normal sample. Default is 1000

### `--max-tmr-cov`
Maximum variant coverage in the tumor sample. Default is 1000

### `--min-fisher`
Minimum phred scaled score for somatic variants. Default is 5

### `--min-str-fisher`
Minimum phred scaled score for STR variants. Default is 25

### `--min-strand-cnt`
Minimum per strand contribution for a variant. Default is 1

### Feature Flags
Use these flags to toggle certain portions of the code to apply different features. By default, these features are off, but by using tags (with no argument following), you can toggle the feature on.

### `--verbose`
Turn on verbose logging for more detailed messages

### `--tenx-mode`
Run Lnacet in 10X Linked Reads mode

### `--active-region-off`
Turn off active region detection

### `--kmer-recovery-on`
Turn on experimental kmer recovery

### `--xa-filter`
Skip reads with XA tag (BWA-mem only)

### `--skip-secondary`
Skip secondary read alignments

### `--extract-pairs`
Extract read pairs for each window

### `--no-contig-check`
Skip checking for same contigs in BAM/CRAMs and reference

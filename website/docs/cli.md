---
id: cli
title: Command Line Reference
---
## pipeline
This is the subcommand that will kick off the tool. This will always follow directly after the call to the Lancet executable

```bash
./Lancet2 pipeline --help
./Lancet2 pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out_prefix
```

## General Options:
### `-h`, `--help`
Bring up the standard help output. This will give examples of options the user can change to customize a run.

### `-v`, `--version`
Print the version information for the build

## Required Arguments:
A standard run of Lancet will provide a tumor bam file, a normal bam file, the reference fasta, and an output path to place the outputted vcf file.
```shell
./lancet pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out.vcf
```

### `-t`, `--tumor`
Provide the path to the tumor BAM/CRAM file(s). Index for this file should also be in the same directory

### `-n`, `--normal`
Provide the path to normal BAM/CRAM file(s). Index for this file should also be in the same directory

### `-r`, `--reference`
Provide the path to the reference fasta file. Index for this file should also be in the same directory

### `-o`, `--out-vcfgz`
Output path to the compressed VCF file

## Optional Arguments:
These arguments allow for more fine-tuned control of the tool. If not provided, default values will be assigned
### `--graphs-dir`
This tag allows you to define the output path for dumping serialized graphs from a run. If this option is not utilized, there will be no outputted graphs.

### Regions
These options will allow you to play around with what the tool looks at.

### `-R`, `--region`
Allows the user to define what region(s) the tool should run on. Should be of format chr:start_pos-end_pos and multiple regions will be separated by a comma.
For example, this will indicate that the tool will run on chromosome 2 between positions 33091000 and 33092000 as well as chromosome 18 between positions 5729000 and 5759000:
```shell
... --region 2:33091000-33092000,18:5729000-5759000 ...
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
Padding to be applied to all input regions. By default, there will be a padding of 500 bp.

### `-w`, `--window-size`
This tag allows you to define the window size used for each microassembly task. By default, the window will be 500 bp.

### `-p`, `--pct-overlap`
Allows the user to define how much overlap there should be between windows. If not specified, the tool will default to 20% overlap

### Parameters
These options allow you to define certain parameters for how the tool performs its variant calling

### `-T`, `--num-threads`
This allows you to define how many async worker threads are used by the tool. If not specified, the tool will default to running with 2 threads.

### `-k`, `--min-kmer`
This allows you to define the minimum length kmers should be for graph nodes. If no length specified, the min kmer length defaults to 31 bp.

### `-K`, `--max-kmer`
This allows you to define the maximum length kmers should be for graph nodes. If no length specified, the max kmer length defaults to 133 bp.

### `-s`, `--kmer-step
This allows you to define the kmer step sizes to try for graph nodes. If not specified, the kmer step size defaults to 4

### `--min-anchor-cov`
Minimum coverage for anchor nodes (source & sink). Default is 5

### `--min-node-cov`
Minimum coverage for all nodes in the graph. Default is 2

### `--max-sample-cov`
Maximum per sample coverage before downsampling. Default is 500

### `--min-alt-qual`
Minimum phred quality supporting ALT allele. Default is 20

### Filters
These options let you apply different filters to apply to the variant caller

### `--min-nml-cov`
Minimum normal coverage. Default is 20

### `--min-tmr-cov`
Minimum tumor coverage. Default is 20

### `--min-odds-ratio`
Minimum VAF odds in tumor vs normal. Default is 10

### `--min-snv-fisher`
Minimum phred scaled fisher score for SNVs

### `--min-indel-fisher`
Minimum phred scaled fisher score for InDels

### `--min-str-fisher`
Minimum phred scaled fisher score for STRs

### Flags
Use these flags to toggle certain portions of the code to apply different features. By default, these features are off, but by using tags (with no argument following), you can toggle the feature on.

### `--verbose`
Turn on verbose logging for more detailed messages

### `--extract-pairs`
Extract all useful read pairs

### `--no-active-region`
Force assemble all windows

### `--no-contig-check`
Skip contig check with reference

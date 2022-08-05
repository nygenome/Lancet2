---
id: cli
title: Command Line Options
---
## pipeline
This is the subcommand that will kick off the tool. This will always follow directly after the call to the Lancet executable
```shell
./lancet pipeline --help
./lancet pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out.vcf
```

## General Options:
### `--help`
Bring up the standard help output. This will give examples of options the user can change to customize a run.

### `--version`
Print the version information for the build

### `--verbose`
Turns on verbose logging and will allow the user to see all log outputs

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

### `-o`, `--out-vcf`
Specify where to place the outputted vcf file.`

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

### 

---
sidebar_position: 1
---
# Basic Usage

Here is a basic run of the Lancet tool from the newly created build directory:

```bash
lancet2 pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/out --region 22 --num-threads 8
```

The command above detects somatic variants in a tumor/normal pair of bam files for chromosome 22 using 8 threads and outputs the VCF file as out.vcf.gz.

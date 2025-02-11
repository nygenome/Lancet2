# Targeted Analysis

Lancet can be run on exome/panel data by providing a preselected list of regions along the genome with a BED file.
The BED file is a tab delimited text file which must contain at least three columns in the following order:
- First column indicates the chromosome of the desired region to analyze
- Second and Third columns are the start and end positions of the region respectively.

`sample.bed`
```bash
chr1   56091000    56092000
chr5   37281200    37291200
chr8   11200000    11300000
```

In the example above the `sample.bed` file can be used to call variants in chromosomes 1, 5, and 8 from
positions 56091000-56092000, 37281200-37291200, and 11200000-11300000 respectively using the following command:

```bash
Lancet2 pipeline --reference ref.fasta \
    --tumor tumor.bam --normal normal.bam \
    --bed-file sample.bed --out-vcfgz output.vcf.gz
```

!!! note "Note"

    Chromosome names in the BED file must match the chromosome names present in the reference FASTA and BAM/CRAM files.
    Though it is not recommened, you can use the `--no-contig-check` flag to force ignore this check.

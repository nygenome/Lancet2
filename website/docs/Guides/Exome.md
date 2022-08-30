# Exome Sequencing

Lancet can be run on a preselected list of sequences from the whole genome by providing a BED file. The BED file is a tab delimited text file which must contain at least three columns in the following order: the first column indicates the chromosome of the desired region to analyze and the second and third columns are the start and end positions of the region respectively.

sample.bed:
```bed
1	56091000	56092000
5	37281200	37291200
8	11200000	11300000	
```

In the example above the ```sample.bed``` file can be used to call variants in chromosomes 1, 5, and 8 from positions 56091000-56092000, 37281200-37291200, and 11200000-11300000 respectively using the following command:

```bash
./lancet2 pipeline -t tumor.bam -n normal.bam -r ref.fasta -o out.vcf -b sample.bed
```
NOTE: Chromosome labels in the BED file must match the chromosome labels present in the genome reference and BAM files.

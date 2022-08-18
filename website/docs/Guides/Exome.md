# Exome Sequencing

Lancet can be run on multiple specific sequences in the full genome by providing a bed file. The bed file should be a tab delimited file with three columns. The first column indicates which chromosome the desired region is in and the second and third columns are the start and end positions of the region respectively.

sample.bed:
```bed
1	56091000	56092000
5	37281200	37291200
8	11200000	11300000	
```

The above example can be given as an argument to the tool and call variants from chromosomes 1, 5, and 8 from positions 56091000-56092000, 37281200-37291200, and 11200000-11300000 respectively.

```bash
./lancet2 pipeline -t tumor.bam -n normal.bam -r ref.fasta -o out.vcf -b sample.bed
```

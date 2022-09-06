# Genome-Wide Scan

For whole-genome sequencing studies it is highly recommended to split the analysis by chromosome and then merge the results. Splitting the work by chromosome will reduce overall runtime and memory requirements to analyze whole-genome data.

```bash
NUMBER_OF_AUTOSOMES=22
for chrom in `seq 1 $NUMBER_OF_AUTOSOMES` X Y; do
 qsub \
 -N lancet_chr${chrom} \
 -cwd \
 -pe smp 8 \
 -q dev.q \
 -j y \
 -b y \
 "lancet2 -t T.bam -n N.bam -r ref.fa --region $chrom --num-threads 8 -o ${chrom}_out"
done

// merge VCF files
```

The previous command shows an example submission of multiple parallel lancet jobs, one for each human chromosome, to the Sun Grid Engine queuing system.

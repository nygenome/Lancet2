# Whole Genome Analysis

For whole-genome sequencing studies it is highly recommended to split the analysis by chromosome and then merge the results.
Splitting the work by chromosome will reduce overall runtime requirements to analyze whole-genome data.

```bash
NUM_CORES=64

for chrom in $(head -24 GRCh38.fasta.fai | cut -f1 | tr '\n' ' ')
do
qsub -N "Lancet2_${chrom}" -cwd -pe smp "${NUM_CORES}" -j y -b y \
"Lancet2 pipeline --threads ${NUM_CORES} \
    --normal normal.bam --tumor tumor.bam \
    --reference GRCh38.fasta --region ${chrom} \
    --out-vcfgz output.${chrom}.vcf.gz"
done

# merge per chromosome VCF files
```

The previous command shows an example submission of multiple parallel lancet jobs,
one for each human chromosome, to the Sun Grid Engine queuing system.

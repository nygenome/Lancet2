# Genome-Wide Scan

Due to its pure local-assembly strategy, Lancet currently has longer runtimes compared to standard alignment-based variant callers. For whole-genome sequencing studies it is highly recommended to split the analysis by chromosome and then merge the results. Splitting the work by chromosome will also reduce the overall memory requirements to analyze the whole-genome data.

```
NUMBER_OF_AUTOSOMES=22
for chrom in `seq 1 $NUMBER_OF_AUTOSOMES` X Y; do
	qsub \
	-N lancet_chr${chrom} \
	-cwd \
	-pe smp 8 \
	-q dev.q \
	-j y \
	-b y \
	"lancet2 -t T.bam -n N.bam -r ref.fa --region $chrom --num-threads 8 -o ${chrom}.vcf"
done

// merge VCF files
```
The previous command shows an exemplary submission of multiple parallel lancet jobs, one for each human chromosome, to the Sun Grid Engine queuing system.

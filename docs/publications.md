---
hide:
  - navigation
---

# Publications

## [Lancet2: Improved and accelerated somatic variant calling with joint multi-sample local assembly graphs](https://www.biorxiv.org/content/10.1101/2025.02.18.638852v2.full)
**Abstract:**

Here, we present Lancet2, an open-source somatic variant caller designed to improve detection of small variants in short-read sequencing data.
Lancet2 introduces significant enhancements, including: 1) Improved variant discovery and genotyping through partial order multiple sequence
alignment of assembled haplotype contigs and re-alignment of sample reads to the best supporting allele. 2) Optimized somatic variant scoring
with Explainable Machine Learning models leading to better somatic filtering throughout the sensitivity scale. 3) Integration with Sequence Tube Map
for enhanced visualization of variants with aligned sample reads in graph space. When benchmarked against enhanced two-tech truth sets generated
using high-coverage short-read (Illumina) and long-read (Oxford Nanopore) data from four well characterized matched tumor/normal cell lines,
Lancet2 outperformed other industry-leading tools in variant calling performance, especially for InDels. In addition, significant runtime
performance improvements compared to Lancet1 (~10x speedup and 50% less peak memory usage) and most other state of the art somatic variant
callers (at least 2x speedup with 8 cores or more) make Lancet2 an ideal tool for accurate and efficient somatic variant calling.

- **Authors** – Rajeeva Lochan Musunuri, Bryan Zhu, Wayne E. Clarke, William Hooper, Timothy Chu, Jennifer Shelton, André Corvelo, Dickson Chung,
                Shreya Sundar, Adam M Novak, Benedict Paten, Michael C. Zody, Nicolas Robine, Giuseppe Narzisi
- **Journal** – bioRxiv, February 2025

## [Somatic variant analysis of linked-reads sequencing data with Lancet](https://academic.oup.com/bioinformatics/article/37/13/1918/5926970)
**Abstract:**

We present a new version of the popular somatic variant caller, Lancet, that supports the analysis of linked-reads sequencing data.
By seamlessly integrating barcodes and haplotype read assignments within the colored De Bruijn graph local-assembly framework,
Lancet computes a barcode-aware coverage and identifies variants that disagree with the local haplotype structure.

- **Authors** – Rajeeva Musunuri, Kanika Arora, André Corvelo, Minita Shah, Jennifer Shelton, Michael C Zody, Giuseppe Narzisi.
- **Journal** – Bioinformatics, Volume 37, Issue 13, July 2021, Pages 1918–1919.

## [Genome-wide somatic variant calling using localized colored de Bruijn graphs](https://www.nature.com/articles/s42003-018-0023-9)
**Abstract:**

Reliable detection of somatic variations is of critical importance in cancer research. Here we present Lancet, an accurate and sensitive
somatic variant caller, which detects SNVs and indels by jointly analyzing reads from tumor and matched normal samples using colored de Bruijn
graphs. We demonstrate, through extensive experimental comparison on synthetic and real whole-genome sequencing datasets, that Lancet has better
accuracy, especially for indel detection, than widely used somatic callers, such as MuTect, MuTect2, LoFreq, Strelka, and Strelka2. Lancet features
a reliable variant scoring system, which is essential for variant prioritization, and detects low-frequency mutations without sacrificing the
sensitivity to call longer insertions and deletions empowered by the local-assembly engine. In addition to genome-wide analysis, Lancet allows
inspection of somatic variants in graph space, which augments the traditional read alignment visualization to help confirm a variant of interest.
Lancet is available as an open-source program at https://github.com/nygenome/lancet.

- **Authors** – Giuseppe Narzisi, André Corvelo, Kanika Arora, Ewa A. Bergmann, Minita Shah, Rajeeva Musunuri, Anne-Katrin Emde, Nicolas Robine, Vladimir Vacic & Michael C. Zody.
- **Journal** – Nature, Communications Biology 1, Article Number: 20, March 2018.

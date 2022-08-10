---
id: publications
title: Publications
---
# Publications

### [Genome-wide somatic variant calling using localized colored de Bruijn graphs](https://www.nature.com/articles/s42003-018-0023-9)

Abstract:

Reliable detection of somatic variations is of critical importance in cancer research. Here we present Lancet, an accurate and sensitive somatic variant caller, which detects SNVs and indels by jointly analyzing reads from tumor and matched normal samples using colored de Bruijn graphs. We demonstrate, through extensive experimental comparison on synthetic and real whole-genome sequencing datasets, that Lancet has better accuracy, especially for indel detection, than widely used somatic callers, such as MuTect, MuTect2, LoFreq, Strelka, and Strelka2. Lancet features a reliable variant scoring system, which is essential for variant prioritization, and detects low-frequency mutations without sacrificing the sensitivity to call longer insertions and deletions empowered by the local-assembly engine. In addition to genome-wide analysis, Lancet allows inspection of somatic variants in graph space, which augments the traditional read alignment visualization to help confirm a variant of interest. Lancet is available as an open-source program at https://github.com/nygenome/lancet.

### [Somatic variant analysis of linked-reads sequencing data with Lancet](https://www.biorxiv.org/content/10.1101/2020.07.04.158063v1)

Abstract:

We present a new version of the popular somatic variant caller, Lancet, that supports the analysis of linked-reads sequencing data. By seamlessly integrating barcodes and haplotype read assignments within the colored De Bruijn graph local-assembly framework, Lancet computes a barcode-aware coverage and identifies variants that disagree with the local haplotype structure.

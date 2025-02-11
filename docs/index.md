---
hide:
  - navigation
---

# Getting Started

## Background

Lancet2 is a command line somatic variant caller (SNVs and InDels) for short
read sequencing data implemented with modern C++. It performs joint multi-sample
localized colored de-bruijn graph assembly for more accurate variant calls,
especially InDels.

## Installation

### Build prerequisites
- [Linux](https://kernel.org/) [x86-64 system](https://en.wikipedia.org/wiki/X86-64)
- [Git](https://command-not-found.com/git), [Make](https://command-not-found.com/make)
- [GCC](https://gcc.gnu.org) (12.x or greater)
- [CMake](https://cmake.org/download) (3.25 or greater)
- [BZip2](https://sourceware.org/bzip2/), [LibLZMA](https://tukaani.org/xz/)

### Build commands
```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc)
```

### Static binary

!!! note "Note"

    It is recommended to build Lancet2 from scratch on the target machine
    where processing is expected to happen for maximum runtime performance.

The simplest way to use `Lancet2` is to download the binary from the [latest available release](https://github.com/nygenome/Lancet2/releases).
The binary from releases is static, with no dependencies and needs only executable permissions before it can be used.

```bash
chmod +x Lancet2
./Lancet2 --help
```

### Docker images

!!! note "Note"

    A CPU that supports the [AVX512 instruction set](https://en.wikipedia.org/wiki/AVX-512) is required
    to use the pre-built public docker images. Custom docker images for older CPUs can be built by the user by
    modifying the `BUILD_ARCH` argument in the [Dockerfile](https://github.com/nygenome/Lancet2/blob/main/Dockerfile).

Public docker images hosted on Google Cloud are available for [recent tagged releases](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet).

## Basic Usage
The following command demonstrates the basic usage of the Lancet2 variant calling pipeline for a tumor and normal bam file pair on chr22.

```bash
Lancet2 pipeline \
    --normal /path/to/normal.bam \
    --tumor /path/to/tumor.bam \
    --reference /path/to/reference.fasta \
    --region "chr22" --num-threads $(nproc) \
    --out-vcfgz /path/to/output.vcf.gz
```

## License

Lancet2 is distributed under the [BSD 3-Clause License](https://github.com/nygenome/Lancet2/blob/main/LICENSE).

## Citing Lancet2

- [Genome-wide somatic variant calling using localized colored de Bruijn graphs](https://www.nature.com/articles/s42003-018-0023-9)
- [Somatic variant analysis of linked-reads sequencing data with Lancet](https://academic.oup.com/bioinformatics/article/37/13/1918/5926970)

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

# Quick Start

## Background

Lancet is a somatic variant caller (SNVs and indels) for short read data. Lancet uses a localized micro-assembly strategy to detect somatic mutation with high sensitivity and accuracy on a tumor/normal pair.
Lancet is based on the colored de Bruijn graph assembly paradigm where tumor and normal reads are jointly analyzed within the same graph. On-the-fly repeat composition analysis and self-tuning k-mer strategy are used together to increase specificity in regions characterized by low complexity sequences. Lancet requires the raw reads to be aligned with BWA (See [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) description for more info). Lancet is implemented in C++.

## Installation

Lancet can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested.

The recommend way to install lancet for most users is:

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build && cmake .. && make
```

You can then optionally add the built binary (`lancet2`) to your `PATH`:

```bash
echo 'export PATH='$(pwd)':$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

To check for a successful install:

```bash
lancet2 --version
```

### Prerequisites to compile from source

* Linux
* x86-64
* CMake (3.25 or greater)
* GCC (12.x or greater)
* BZip2
* LibLZMA
* Git
* Make

### Docker

Pre-built docker images for Lancet2 are available [here](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet). Note that the public docker image requires [Cascade Lake CPU](https://en.wikichip.org/wiki/intel/microarchitectures/cascade_lake) or newer

You can also build a new image from the Dockerfile provided in the repository.

## Basic Usage

Here is a basic run of the Lancet tool from the newly created build directory:

```bash
Lancet2 pipeline \
    -t /path/to/tumor.bam \
    -n /path/to/normal.bam \
    -r /path/to/ref.fasta \
    --region 22 --num-threads 8 \
    -o /path/to/out
```

The command above detects somatic variants in a tumor/normal pair of bam files for chromosome 22 using 8 threads and outputs the VCF file as out.vcf.gz.

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

In addition to variant calling accuracy and improved somatic filtering, Lancet2 has significant runtime performance improvements compared to Lancet1 (upto ∼10x speedup and 50% less peak memory usage)

## Installation

### Build prerequisites
- [Linux](https://kernel.org/) or [macOS](https://www.apple.com/macos/) (x86-64 or ARM64 architectures)
- [Git](https://command-not-found.com/git), [Make](https://command-not-found.com/make)
- [GCC](https://gcc.gnu.org) (12.x or greater) or [Clang](https://clang.llvm.org) (14.x or greater)
- [CMake](https://cmake.org/download) (3.25 or greater)
- [BZip2](https://sourceware.org/bzip2/), [LibLZMA](https://tukaani.org/xz/)
- [CURL](https://curl.se/) and [OpenSSL](https://www.openssl.org/) (optional, required for native Cloud I/O)

### Build commands
```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc)
```

!!! note "macOS Homebrew Configuration"

    To fix CMake failing to find Homebrew-installed BZip2 on macOS, explicitly set the `CMAKE_PREFIX_PATH` to the Homebrew prefix when running CMake. Homebrew often skips linking `bzip2` into the shared prefix to prevent conflicts with system libraries.
    
    For Apple Silicon (M1/M2/M3):
    ```
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/bzip2 ..
    ```
    
    For Intel Mac:
    ```
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr/local/opt/bzip2 ..
    ```

### Cloud I/O Support (GCS, S3, HTTP/S, FTP/S)
Native support for network streaming from Google Cloud Storage (`gs://`), Amazon S3 (`s3://`), and standard web endpoints (`http(s)://` & `ftp(s)://`) can be enabled by setting `-DLANCET_ENABLE_CLOUD_IO=ON` during CMake configuration. This feature is opt-in and is only permitted when configuring a dynamically linked build (`-DLANCET_BUILD_STATIC=OFF`). It requires the system to have `libcurl` and `openssl` installed. For exact usage instructions, required environment variables, and authentication configuration, please refer to the [Native Cloud Streaming Guide](guides/cloud_streaming.md).
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON ..
make -j$(nproc)
```

### Static binary

!!! note "Note"

    It is recommended to build Lancet2 from scratch on the target machine
    where processing is expected to happen for maximum runtime performance.

If you have a Linux based operating system and a [CPU that supports AVX2 instructions](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX2).
The simplest way to use `Lancet2` is to download the binary from the [latest available release](https://github.com/nygenome/Lancet2/releases).
The binary from releases is static, with no dependencies and needs only executable permissions before it can be used. 

```bash
chmod +x Lancet2
./Lancet2 --help
```

> [!WARNING]
> Because Cloud I/O fundamentally requires dynamically linking the underlying host OS's network stack (`libcurl` / `openssl`), these statically compiled release binaries **do not** support Native Cloud Streaming (`gs://`, `s3://`). If you require direct cloud network capabilities without compiling from source, please use the official **Docker** containers (Conda packaging is actively in-progress), which organically provide dynamic networking dependencies out-of-the-box.

### Docker images

!!! note "Note"

    A CPU that supports the [AVX2 instruction set](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX2) is required
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

See the [Scoring Somatic Variants](guides/scoring_somatic_variants.md) guide for more information on how
to score and filter somatic variants using explainable machine learning models.

## License

Lancet2 is distributed under the [BSD 3-Clause License](https://github.com/nygenome/Lancet2/blob/main/LICENSE).

## Citing Lancet2

- [Lancet2: Improved and accelerated somatic variant calling with joint multi-sample local assembly graphs](https://www.biorxiv.org/content/10.1101/2025.02.18.638852v2.full)
- [Somatic variant analysis of linked-reads sequencing data with Lancet](https://academic.oup.com/bioinformatics/article/37/13/1918/5926970)
- [Genome-wide somatic variant calling using localized colored de Bruijn graphs](https://www.nature.com/articles/s42003-018-0023-9)

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

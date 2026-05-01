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

### Pre-built packages (Recommended)
Lancet2 packages with **full Cloud I/O support** (`s3://`, `gs://`, `http(s)://`, `ftp(s)://`) are published to [prefix.dev/channels/lancet2](https://prefix.dev/channels/lancet2). Install using your preferred package manager:

| Package Manager | Install Command |
|---|---|
| [Pixi](https://pixi.sh/) (recommended) | `pixi global install --channel https://prefix.dev/channels/lancet2 lancet2` |
| [Conda](https://docs.conda.io/) | `conda install -c https://prefix.dev/channels/lancet2 lancet2` |
| [Mamba](https://mamba.readthedocs.io/) | `mamba install -c https://prefix.dev/channels/lancet2 lancet2` |

Development builds are published automatically on every commit. To install a specific stable release, pin the version:
```bash
pixi global install --channel https://prefix.dev/channels/lancet2 'lancet2==v2.9.0'
```

### Docker images

!!! note "Note"

    A CPU that supports the [AVX2 instruction set](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX2) is required
    to use the pre-built public docker images. Custom docker images for older CPUs can be built by the user by
    modifying the `BUILD_ARCH` argument in the [Dockerfile](https://github.com/nygenome/Lancet2/blob/main/Dockerfile).

Public docker images hosted on Google Cloud are available for [recent tagged releases](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet).

### Build from source

!!! note "Note"

    Building from source on the target machine is recommended for maximum runtime performance.

!!! warning "Cloud support in static builds"

    Static builds (the default) **do not** support Cloud Streaming (`gs://`, `s3://`, `http(s)://`, `ftp(s)://`) because cloud I/O requires dynamic linking of the host OS's network stack (`libcurl` / `openssl`).
    Use the pre-built packages or Docker images instead, or build with `-DLANCET_ENABLE_CLOUD_IO=ON` (see below).

#### Build prerequisites
- [Linux](https://kernel.org/) or [macOS](https://www.apple.com/macos/) (x86-64 or ARM64 architectures)
- [Make](https://command-not-found.com/make)
- [GCC](https://gcc.gnu.org) (12.x or greater) or [Clang](https://clang.llvm.org) (14.x or greater)
- [CMake](https://cmake.org/download) (3.25 or greater)
- [BZip2](https://sourceware.org/bzip2/), [LibLZMA](https://tukaani.org/xz/)
- [CURL](https://curl.se/) and [OpenSSL](https://www.openssl.org/) (optional, required for Cloud I/O)

#### Build commands
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

#### Cloud I/O Support (GCS, S3, HTTP/S, FTP/S)
Cloud streaming can be enabled by setting `-DLANCET_ENABLE_CLOUD_IO=ON` during CMake configuration. This is opt-in and requires dynamic linkage (`-DLANCET_BUILD_STATIC=OFF`) with `libcurl` and `openssl` installed. For usage details, see the [Cloud Streaming Guide](guides/cloud_streaming.md).
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON ..
make -j$(nproc)
```

## Basic Usage

### Tumor-Normal Somatic Calling

The primary supported workflow. Variants are classified as `CASE` (tumor-only), `CTRL` (normal-only), or `SHARED` (both) in the VCF INFO field.

```bash
Lancet2 pipeline \
    --normal /path/to/normal.bam \
    --tumor /path/to/tumor.bam \
    --reference /path/to/reference.fasta \
    --region "chr22" --num-threads $(nproc) \
    --out-vcfgz /path/to/output.vcf.gz
```

!!! note "Why `CASE`/`CTRL` instead of `TUMOR`/`NORMAL`?"
    Lancet2 uses **case/control** terminology to generalize beyond tumor-normal analysis (e.g., treated vs. untreated, responder vs. non-responder). The `--normal` and `--tumor` CLI flags are **permanent and first-class** — tumor-normal somatic calling is the primary supported workflow. See the [CLI Reference](reference.md#datasets) for details.

### Single-Sample Mode

Omit `--tumor` to run on a single sample. Lancet2 generates raw variant candidates across the full allele spectrum (germline, mosaic, artifact) without somatic state classification — no `SHARED`/`CTRL`/`CASE` tags appear in the VCF.

```bash
Lancet2 pipeline \
    --normal /path/to/sample.bam \
    --reference /path/to/reference.fasta \
    --region "chr22" --num-threads $(nproc) \
    --out-vcfgz /path/to/output.vcf.gz
```

### Multi-Sample Mode

Additional samples can be added using [`-s,--sample`](reference.md#datasets) alongside the standard flags. See [Multi-Sample & Germline Mode](guides/architecture.md#multi-sample-germline-mode) for details.

!!! warning "Experimental — no pre-trained ML models"
    Single-sample and multi-sample modes produce raw variant candidates. No pre-trained ML models are currently provided for filtering in these modes — variant calls require custom downstream filtering. The only pre-trained model available is for the standard tumor-normal somatic workflow (v2.8.7 compatible). See [Scoring Somatic Variants](guides/scoring_somatic_variants.md) for details.

## License

Lancet2 is distributed under the [BSD 3-Clause License](https://github.com/nygenome/Lancet2/blob/main/LICENSE).

## Citing Lancet2

- [Lancet2: Improved and accelerated somatic variant calling with joint multi-sample local assembly graphs](https://academic.oup.com/nargab/article/8/2/lqag036/8625920)
- [Somatic variant analysis of linked-reads sequencing data with Lancet](https://academic.oup.com/bioinformatics/article/37/13/1918/5926970)
- [Genome-wide somatic variant calling using localized colored de Bruijn graphs](https://www.nature.com/articles/s42003-018-0023-9)

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

---
hide:
  - navigation
---

# Getting Started

## Background

Lancet2 is an accurate variant caller (SNVs and InDels) for short read
sequencing data, supporting flexible single-sample through multi-sample
analysis via localized colored de Bruijn graph assembly.

In addition to variant calling accuracy and improved somatic filtering, Lancet2 has significant
runtime performance improvements compared to Lancet1 (up to ~10× speedup and ~50% less peak memory usage).

## Installation

Install a pre-built package using [Pixi](https://pixi.sh/), [Conda](https://docs.conda.io/), or [Mamba](https://mamba.readthedocs.io/).
Packages with **Cloud I/O support** (`s3://`, `gs://`, `http(s)://`, `ftp(s)://`) are published to [prefix.dev/channels/lancet2](https://prefix.dev/channels/lancet2):

| Package Manager | Install Command |
|---|---|
| [Pixi](https://pixi.sh/) (recommended) | `pixi global install --channel https://prefix.dev/channels/lancet2 lancet2` |
| [Conda](https://docs.conda.io/) | `conda install -c https://prefix.dev/channels/lancet2 lancet2` |
| [Mamba](https://mamba.readthedocs.io/) | `mamba install -c https://prefix.dev/channels/lancet2 lancet2` |

```bash
# Find the latest stable version at https://prefix.dev/channels/lancet2/packages/lancet2
# Stable releases: clean numbers (e.g., 2.9.0)
# Dev builds: version_branch_hash (e.g., 2.9.0_main_441d081927)
export LANCET_VERSION="2.9.0"
pixi global install --channel https://prefix.dev/channels/lancet2 "lancet2==${LANCET_VERSION}"
```

!!! note "Hardware and platform support"
    Pre-built packages are available for linux-64 (x86-64-v3, Haswell 2013+)
    and macOS ARM64 (Apple Silicon M1+). Intel Macs can
    [build from source](guides/installation.md#standalone-build).
    See the [Installation Guide](guides/installation.md) for Docker,
    build-from-source, and Cloud I/O instructions.

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

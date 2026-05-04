# Lancet2

Lancet2 is an accurate variant caller (SNVs and InDels) for short read
sequencing data, supporting flexible single-sample through multi-sample
analysis via localized colored de Bruijn graph assembly.

In addition to variant calling accuracy and improved somatic filtering, Lancet2 has significant
runtime performance improvements compared to Lancet1 (up to ~10× speedup and ~50% less peak memory usage).

[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?mLabel=Documentation&style=flat)](https://nygenome.github.io/Lancet2)
[![GitHub Release](https://img.shields.io/github/v/release/nygenome/Lancet2?include_prereleases&sort=semver&display_name=release&style=flat)](https://github.com/nygenome/Lancet2/releases)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Installation

### Pre-built packages (Recommended)
Lancet2 packages with **Cloud I/O support** (`s3://`, `gs://`, `http(s)://`, `ftp(s)://`) are published to [prefix.dev/channels/lancet2](https://prefix.dev/channels/lancet2).
Install using your preferred package manager:

| Package Manager | Install Command |
|---|---|
| [Pixi](https://pixi.sh/) (recommended) | `pixi global install --channel https://prefix.dev/channels/lancet2 lancet2` |
| [Conda](https://docs.conda.io/) | `conda install -c https://prefix.dev/channels/lancet2 lancet2` |
| [Mamba](https://mamba.readthedocs.io/) | `mamba install -c https://prefix.dev/channels/lancet2 lancet2` |

Development builds are published on every push to `main`.
Stable releases are published with every tagged release.
To install a specific stable release, pin the version:

```bash
# Stable releases: clean numbers (e.g., 2.9.0)
# Dev builds: version_branch_hash (e.g., 2.9.0_main_441d081927)
# Find versions at https://prefix.dev/channels/lancet2/packages/lancet2
export LANCET_VERSION="2.9.0"
pixi global install --channel https://prefix.dev/channels/lancet2 "lancet2==${LANCET_VERSION}"
```

Pre-built packages are available for **linux-64** (x86-64-v3, Haswell 2013+) and **macOS ARM64** (Apple Silicon M1+).

### Docker images

Pre-built Docker images are available from [Google Artifact Registry](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet).
A CPU with x86-64-v3 support (Haswell 2013+: AVX2, BMI2, FMA) is required.
Docker tags follow the same convention: stable releases use clean version numbers (`X.Y.Z`),
dev builds include branch and commit hash (`X.Y.Z_main_<hash>`).

```bash
docker pull us-central1-docker.pkg.dev/nygc-app-c-148c/lancet-public/lancet:X.Y.Z
```

### Build from source

See the [Installation Guide](https://nygenome.github.io/Lancet2/latest/guides/installation/#build-from-source) for build-from-source instructions, including static and Cloud I/O build variants.

## Quick Start

Tumor-normal somatic calling — the primary supported workflow. In the VCF,
`--tumor` maps to `CASE` and `--normal` maps to `CTRL`
([why?](https://nygenome.github.io/Lancet2/latest/#basic-usage)):

```bash
Lancet2 pipeline \
    --normal /path/to/normal.bam \
    --tumor /path/to/tumor.bam \
    --reference /path/to/reference.fasta \
    --region "chr22" --num-threads $(nproc) \
    --out-vcfgz /path/to/output.vcf.gz
```

Lancet2 also supports **single-sample** (germline/mosaic) and **multi-sample**
modes. See the [full documentation](https://nygenome.github.io/Lancet2/) for
details.

## Documentation

See [full documentation](https://nygenome.github.io/Lancet2/) for
[VCF output format](https://nygenome.github.io/Lancet2/guides/vcf_output/),
[CLI reference](https://nygenome.github.io/Lancet2/reference/), and
[pipeline architecture](https://nygenome.github.io/Lancet2/guides/architecture/).

## Citing

See [publications](https://nygenome.github.io/Lancet2/publications) associated with Lancet.

## License

Lancet2 is distributed under the [BSD 3-Clause License](LICENSE).

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

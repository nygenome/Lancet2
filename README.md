# Lancet2

Lancet2 is a command line somatic variant caller (SNVs and InDels) for short
read sequencing data implemented with modern C++. It performs joint multi-sample
localized colored de-bruijn graph assembly for more accurate variant calls,
especially InDels.

[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?mLabel=Documentation&style=flat)](https://nygenome.github.io/Lancet2)
[![GitHub Release](https://img.shields.io/github/v/release/nygenome/Lancet2?include_prereleases&sort=semver&display_name=release&style=flat)](https://github.com/nygenome/Lancet2/releases)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

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

## Documentation

Documentation for Lancet2 is hosted on [GitHub pages](https://nygenome.github.io/Lancet2/).

## Citing

See [publications](https://nygenome.github.io/Lancet2/publications) associated with Lancet.

## License

Lancet2 is distributed under the [BSD 3-Clause License](LICENSE).

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

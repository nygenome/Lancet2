# Lancet2

> Micro-assembly based somatic variant caller

***Note: Lancet2 is currently under development and is not yet ready for production use.***

[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?mLabel=Documentation&style=flat)](https://nygenome.github.io/Lancet2)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Installation

### Pre-requisites to compile from source

* [Linux](https://kernel.org/) [x86-64 system](https://en.wikipedia.org/wiki/X86-64)
* [Autoconf](https://command-not-found.com/autoconf)
* [Git](https://command-not-found.com/git)
* [Make](https://command-not-found.com/make)
* [Ninja](https://command-not-found.com/ninja-build)
* [CMake](https://cmake.org/download) (3.25 or greater)
* [GCC](https://gcc.gnu.org) (12.x or greater)

### Build commands

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -GNinja .. && ninja -v
```

## Documentation

Documentation for Lancet2 is hosted on [GitHub pages](https://nygenome.github.io/Lancet2/).

## Authors

Rajeeva Musunuri, Bryan Zhu and Giuseppe Narzisi

## Citing

See [publications](https://nygenome.github.io/Lancet2/docs/publications) associated with Lancet.

## License

Lancet2 is distributed under the [BSD 3-Clause License](LICENSE).

## Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov)) under the NCI U01
award [1U01CA253405-01A1](https://reporter.nih.gov/project-details/10304730).

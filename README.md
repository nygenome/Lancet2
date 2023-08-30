# Lancet2

> Micro-assembly based somatic variant caller

***Note: Lancet2 is currently under development and is not yet ready for production use.***

[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?mLabel=Documentation&style=flat)](https://nygenome.github.io/Lancet2)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Docker builds
* [Docker images](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet) are available for every tagged release.
* Please note that the public docker image requires [Cascade Lake CPU](https://en.wikichip.org/wiki/intel/microarchitectures/cascade_lake) or newer.
* If needed, users can build an image for any 64-bit Linux machine using the [Dockerfile](Dockerfile)

## Installation

### Pre-requisites to compile from source

* [Linux](https://kernel.org/)
* [x86-64](https://en.wikipedia.org/wiki/X86-64)
* [BZip2](https://sourceware.org/bzip2/)
* [LibLZMA](https://tukaani.org/xz/)
* [Git](https://command-not-found.com/git)
* [Make](https://command-not-found.com/make)
* [CMake](https://cmake.org/download) (3.25 or greater)
* [GCC](https://gcc.gnu.org) (12.x or greater)

### Build commands

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc)
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

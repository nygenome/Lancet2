# Lancet2

> Microassembly based somatic variant caller

[![DockerBuildAndPush](https://github.com/nygenome/Lancet2/actions/workflows/main.yml/badge.svg)](https://github.com/nygenome/Lancet2/pkgs/container/lancet2)
[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?label=Documentation&style=flat)](https://nygenome.github.io/Lancet2)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## **NOTE**

Lancet2 is currently under development and is not ready for production use.

## Installation

#### Dependencies

* [CMake](https://cmake.org/download/) >= 3.14.x
* [C++ compiler](https://en.cppreference.com/w/cpp/compiler_support#C.2B.2B17_features) with support
  for [ISO C++17 standard](https://en.cppreference.com/w/cpp/17)
    - Tested with GCC >= 8.x, Clang >= 7.x
* [zlib](https://github.com/madler/zlib), [bzip2](https://github.com/enthought/bzip2-1.0.6)
  , [liblzma](https://tukaani.org/xz/), [cURL](https://curl.haxx.se) and [OpenSSL](https://www.openssl.org)

#### Build commands

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake .. && make
```

## Documentation

Documentation for Lancet2 is hosted on [GitHub pages](https://nygenome.github.io/Lancet2/).

## Support

Please report any bugs or feature requests to the [Lancet2 issue tracker](https://github.com/nygenome/Lancet2/issues).

## Authors

Rajeeva Musunuri, Bryan Zhu and Giuseppe Narzisi

## Citing

See [publications](https://nygenome.github.io/Lancet2/docs/publications) associated with Lancet.

## License

Lancet2 is distributed under the [BSD 3-Clause License](LICENSE).

## Funding
Informatics Technology for Cancer Research (<a href="https://itcr.cancer.gov">ITCR</a>) under the NCI U01 award <a href="https://reporter.nih.gov/project-details/10304730">1U01CA253405-01A1</a>.

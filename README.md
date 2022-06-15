# Lancet v2.x
> Microassembly based somatic variant caller

[![Build Status](https://img.shields.io/travis/com/omicsnut/v2_lancet/master.svg?label=Linux/MacOS&style=flat)](https://travis-ci.com/omicsnut/v2_lancet)
[![Codecov](https://codecov.io/gh/omicsnut/v2_lancet/branch/master/graph/badge.svg)](https://codecov.io/gh/omicsnut/v2_lancet)
[![Documentation](https://img.shields.io/badge/Documentation-latest-blue.svg?label=API%20docs&style=flat)](https://nygenome.github.io/Lancet2)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Installation
#### Dependencies
* [CMake](https://cmake.org/download/) >= 3.14.x
* [C++ compiler](https://en.cppreference.com/w/cpp/compiler_support#C.2B.2B17_features) with support for [ISO C++17 standard](https://en.cppreference.com/w/cpp/17)
    - Tested with GCC >= 8.x, Clang >= 7.x
* [zlib](https://github.com/madler/zlib), [bzip2](https://github.com/enthought/bzip2-1.0.6), [liblzma](https://tukaani.org/xz/), [cURL](https://curl.haxx.se) and [OpenSSL](https://www.openssl.org)

#### Build commands
```bash
git clone https://github.com/omicsnut/v2_lancet
cd lancet && mkdir build && cd build
cmake .. && make
```

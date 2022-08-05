---
sidebar_position: 1
---

# Getting Started

Install the Lancet repository to your machine:

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build
cmake .. && make
```

Here is a basic run of the Lancet tool from the newly created build directory:

```bash
./lancet pipeline -t /path/to/tumor.bam -n /path/to/normal.bam -r /path/to/ref.fasta -o /path/to/output.vcf
```

## Prerequisites:
* CMake >= 3.14.x
* C++ compiler with support for ISO C++17 standard
* zlib, bzip2, liblzma, cURL and OpenSSL


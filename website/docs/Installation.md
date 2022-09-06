# Installation

Lancet can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested.

The recommend way to install lancet for most users is:

```bash
git clone https://github.com/nygenome/Lancet2.git
cd Lancet2 && mkdir build && cd build && cmake .. && make
```

You can then optionally add the built binary (`lancet2`) to your `PATH`:

```bash
echo 'export PATH='$(pwd)':$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

To check for a successful install:

```bash
lancet2 --version
```

## Requirements & Dependencies

* CMake >= 3.14.x
* C++ compiler with support for ISO C++17 standard
* zlib, bzip2, liblzma, cURL and OpenSSL (Needed to build htslib)

## Docker

Pre-build docker images for Lancet2 are available on [DockerHub](https://hub.docker.com/r/rmusunuri/lancet2)

```bash
docker pull rmusunuri/lancet2:gcc
docker run rmusunuri/lancet2:gcc lancet2 --version
```

You can also build a new image from the Dockerfile provided in the repository:

```bash
git clone https://github.com/nygenome/Lancet2.git && cd Lancet2
docker build --platform linux/amd64 -t lancet2 -f docker/gcc.Dockerfile .
docker run --platform linux/amd64 -it lancet2 lancet2 --version
```

## Singularity/Apptainer

To build a [Singularity/Apptainer](https://apptainer.org/) container directly from the DockerHub images use the following command

```bash
singularity build lancet2.sif docker://rmusunuri/lancet2:gcc
```

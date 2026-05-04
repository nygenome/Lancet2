# Installation

## Pre-built Packages (Recommended)

Lancet2 packages with **Cloud I/O support** (`s3://`, `gs://`, `http(s)://`, `ftp(s)://`) are published to [prefix.dev/channels/lancet2](https://prefix.dev/channels/lancet2).
Install using your preferred package manager:

| Package Manager | Install Command |
|---|---|
| [Pixi](https://pixi.sh/) (recommended) | `pixi global install --channel https://prefix.dev/channels/lancet2 lancet2` |
| [Conda](https://docs.conda.io/) | `conda install -c https://prefix.dev/channels/lancet2 lancet2` |
| [Mamba](https://mamba.readthedocs.io/) | `mamba install -c https://prefix.dev/channels/lancet2 lancet2` |

### Choosing a version

Stable releases use clean version numbers (e.g., `2.9.0`).
Development builds include the branch name and commit hash (e.g., `2.9.0_main_441d081927`).
Browse available versions at [prefix.dev/channels/lancet2](https://prefix.dev/channels/lancet2/packages/lancet2).

```bash
# Find the latest stable version at https://prefix.dev/channels/lancet2/packages/lancet2
# Stable releases: clean numbers (2.9.0)
# Dev builds: version_branch_hash (2.9.0_main_441d081927)
export LANCET_VERSION="2.9.0"
pixi global install --channel https://prefix.dev/channels/lancet2 "lancet2==${LANCET_VERSION}"
```

!!! note "Development builds"
    Dev builds are published on every push to `main`. They reflect the
    latest state of development and may contain untested changes. For
    production use, always pin a stable release version.

### Hardware requirements

| Platform | Baseline | Pre-built packages | Build from source |
|:---------|:---------|:-------------------|:------------------|
| linux-64 | x86-64-v3 (Haswell 2013+: AVX2, BMI2, FMA, LZCNT) | ✅ | ✅ |
| osx-arm64 | Apple Silicon (M1+) | ✅ | ✅ |
| osx-64 | Intel Mac (Haswell 2013+) | ❌ No packages built | ✅ Untested but should work |

Pre-Haswell CPUs will crash with `SIGILL` on the first AVX2 instruction.
To target different hardware, [build from source](#build-from-source)
with `-DLANCET_NATIVE_BUILD=ON` or edit `cmake/compiler_flags.cmake`.

---

## Docker

Pre-built Docker images are available from [Google Artifact Registry](https://console.cloud.google.com/artifacts/docker/nygc-app-c-148c/us-central1/lancet-public/lancet).
Same x86-64-v3 hardware requirement as the pre-built packages above.

Docker image tags follow the same convention as the prefix.dev packages:

- **Stable releases** use clean version numbers: `lancet:X.Y.Z` (e.g., `lancet:2.9.0`)
- **Dev builds** include the branch and commit hash: `lancet:X.Y.Z_main_<hash>` (e.g., `lancet:2.9.0_main_441d081927`)

```bash
# Pull a stable release (replace X.Y.Z with the desired version)
docker pull us-central1-docker.pkg.dev/nygc-app-c-148c/lancet-public/lancet:X.Y.Z
```

```bash
# Run (replace X.Y.Z with your installed version)
docker run --rm -v /path/to/data:/data \
  us-central1-docker.pkg.dev/nygc-app-c-148c/lancet-public/lancet:X.Y.Z \
  Lancet2 pipeline \
    --normal /data/normal.bam --tumor /data/tumor.bam \
    --reference /data/ref.fasta --region "chr22" \
    --out-vcfgz /data/output.vcf.gz
```

The [Dockerfile](https://github.com/nygenome/Lancet2/blob/main/Dockerfile)
can be customized for different architectures by setting
`-DLANCET_NATIVE_BUILD=ON` in the CMake configure step.

---

## Build from Source

### Using pixi (Recommended — Linux only)

The [pixi](https://pixi.sh) development environment provides the exact
toolchain used in CI (Clang, Ninja, CMake, and all dependencies).

```bash
git clone https://github.com/nygenome/Lancet2.git && cd Lancet2
curl -fsSL https://pixi.sh/install.sh | bash
```

**Static binary (default, no Cloud I/O):**

```bash
pixi run cmake -GNinja -B build \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DLANCET_BUILD_STATIC=ON -DLANCET_ENABLE_CLOUD_IO=OFF
pixi run cmake --build build
# Binary: build/Lancet2
```

**Dynamic binary with Cloud I/O (S3, GCS, HTTP/FTP):**

```bash
pixi run cmake -GNinja -B build-cloud \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON
pixi run cmake --build build-cloud
# Binary: build-cloud/Lancet2
```

!!! warning "Linux only"
    The pixi environment is configured for `linux-64` only. macOS users
    should use the [standalone build](#standalone-build) below.

### Standalone Build

For macOS or Linux without pixi. Requires:

- GCC ≥ 12 or Clang ≥ 14 (Clang recommended)
- CMake ≥ 3.25
- Ninja (recommended) or Make
- zlib, bzip2, liblzma
- libcurl, openssl (only for Cloud I/O)

**Static binary (default, no Cloud I/O):**

```bash
git clone https://github.com/nygenome/Lancet2.git && cd Lancet2
mkdir build && cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release \
  -DLANCET_BUILD_STATIC=ON -DLANCET_ENABLE_CLOUD_IO=OFF ..
ninja
```

**Dynamic binary with Cloud I/O:**

```bash
git clone https://github.com/nygenome/Lancet2.git && cd Lancet2
mkdir build-cloud && cd build-cloud
cmake -GNinja -DCMAKE_BUILD_TYPE=Release \
  -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON ..
ninja
```

Or with Make if Ninja is not available:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

!!! note "Static vs Dynamic builds"
    Static linking (`LANCET_BUILD_STATIC=ON`, the default) produces a
    single portable binary but cannot include Cloud I/O support because
    cloud streaming requires dynamic linking of `libcurl` and `openssl`.
    Use the pre-built packages if you need Cloud I/O without building
    from source.

See [Native Cloud Streaming](cloud_streaming.md) for authentication
setup and usage details.

### macOS Notes

!!! tip "Homebrew BZip2"
    If CMake fails to find BZip2 on macOS with Homebrew, set:
    ```bash
    export BZIP2_ROOT=$(brew --prefix bzip2)
    ```
    Then re-run the CMake configure step.

macOS does not support static builds (`LANCET_BUILD_STATIC` is forced
`OFF` on Darwin). All macOS builds are dynamically linked.

!!! note "Intel Mac support"
    Building from source on Intel Macs (Haswell 2013+) should work but
    is not tested in CI. Pre-built packages are only available for
    Apple Silicon (osx-arm64).

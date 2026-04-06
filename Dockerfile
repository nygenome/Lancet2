FROM gcc:15 AS builder
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get upgrade --yes --no-install-recommends && \
    apt-get install --yes --no-install-recommends ca-certificates libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev git make unzip wget && \
    wget -cq "https://github.com/Kitware/CMake/releases/download/v4.3.1/cmake-4.3.1-linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.13.1/ninja-linux.zip" && \
    /bin/bash "cmake-4.3.1-linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    unzip -d /usr/bin ninja-linux.zip && rm -f cmake-4.3.1-linux-x86_64.sh ninja-linux.zip

ARG BUILD_ARCH="x86-64-v3"
COPY . /Lancet2
RUN cd /Lancet2 && mkdir build && cd build && \
    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_ARCH=${BUILD_ARCH} -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON .. && ninja -v

# Note: We must use a Debian-based slim image instead of Alpine here.
# Cloud I/O natively links dynamically against libcurl and openSSL, which natively invoke glibc hooks.
# Trying to execute this binary on Alpine's lightweight musl libc will result in terminal execution crashes.
FROM debian:trixie-slim
RUN DEBIAN_FRONTEND="noninteractive" apt-get update && \
    apt-get install --yes --no-install-recommends ca-certificates bash libbz2-dev liblzma-dev curl libcurl4 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
COPY --from=builder /Lancet2/build/Lancet2 /usr/bin/Lancet2

FROM gcc:14 AS builder
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get upgrade --yes --no-install-recommends && \
    apt-get install --yes --no-install-recommends ca-certificates libbz2-dev liblzma-dev git make unzip wget && \
    wget -cq "https://github.com/Kitware/CMake/releases/download/v3.29.5/cmake-3.29.5-linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.12.1/ninja-linux.zip" && \
    /bin/bash "cmake-3.29.5-linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    unzip -d /usr/bin ninja-linux.zip && rm -f cmake-3.29.5-linux-x86_64.sh ninja-linux.zip

ARG BUILD_ARCH="cascadelake"
RUN DEBIAN_FRONTEND="noninteractive" git clone https://github.com/nygenome/Lancet2.git && cd Lancet2 && mkdir build && \
    cd build && cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_ARCH=${BUILD_ARCH} .. && ninja -v

FROM alpine:latest
RUN apk update && apk add --no-cache bash curl
COPY --from=builder /Lancet2/build/Lancet2 /usr/bin/Lancet2

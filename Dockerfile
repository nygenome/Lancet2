FROM gcc:13 AS builder
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get upgrade --yes --no-install-recommends && \
    apt-get install --yes --no-install-recommends ca-certificates unzip git make wget gcc g++ && \
    wget -cq "https://github.com/Kitware/CMake/releases/download/v3.26.4/cmake-3.26.4-linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.11.1/ninja-linux.zip" && \
    /bin/bash "cmake-3.26.4-linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    unzip -d /usr/bin ninja-linux.zip && rm -f cmake-3.26.4-linux-x86_64.sh ninja-linux.zip

RUN DEBIAN_FRONTEND="noninteractive" git clone https://github.com/nygenome/Lancet2.git && cd /Lancet2 && mkdir build && \
    cd build && cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_ARCH="haswell" .. && ninja -v

FROM alpine:latest
RUN apk update && apk add --no-cache bash curl
COPY --from=builder /Lancet2/build/Lancet2 /usr/bin/Lancet2

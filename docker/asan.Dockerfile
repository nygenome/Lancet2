FROM clangbuiltlinux/ubuntu:llvm12-latest

LABEL org.opencontainers.image.authors="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN apt-get update && apt-get -y --no-install-recommends upgrade && \
    apt-get install -y --no-install-recommends zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev && \
    apt-get install -y --no-install-recommends unzip git make wget && \
    wget -cq "https://github.com/Kitware/CMake/releases/download/v3.23.2/cmake-3.23.2-linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.11.0/ninja-linux.zip" && \
    /bin/bash "cmake-3.23.2-linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    unzip -d /usr/bin ninja-linux.zip && apt-get clean && apt-get purge && \
    rm -rf "/var/lib/apt/lists/*" "/tmp/*" "/var/tmp/*" "cmake-3.23.2-linux-x86_64.sh" "ninja-linux.zip"

ENV CC  clang-12
ENV CXX clang++-12

RUN mkdir -p /usr/src && git clone https://github.com/nygenome/Lancet2.git /usr/src/lancet2 && \
    mkdir -p /usr/src/lancet2/build && cd /usr/src/lancet2/build && \
    cmake -GNinja -DCMAKE_C_COMPILER=clang-12 -DCMAKE_CXX_COMPILER=clang++-12 \
          -DCMAKE_BUILD_TYPE=Debug -DLANCET_SANITIZER=Address .. && ninja -v && cp lancet2 /usr/bin/

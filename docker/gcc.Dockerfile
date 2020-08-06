FROM gcc:10

MAINTAINER Rajeeva Musunuri <rmusunuri@nygenome.org>

RUN apt-get update && apt-get -y --no-install-recommends upgrade && \
    apt-get install -y --no-install-recommends zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev && \
    apt-get install -y --no-install-recommends unzip git make && \
    wget -cq "https://github.com/Kitware/CMake/releases/download/v3.18.1/cmake-3.18.1-Linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.10.0/ninja-linux.zip" && \
    /bin/bash "cmake-3.18.1-Linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    unzip -d /usr/bin ninja-linux.zip && apt-get clean && apt-get purge && \
    rm -rf "/var/lib/apt/lists/*" "/tmp/*" "/var/tmp/*" "cmake-3.18.1-Linux-x86_64.sh" "ninja-linux.zip"

ENV CC  gcc
ENV CXX g++

RUN mkdir -p /usr/src && git clone https://github.com/omicsnut/v2_lancet.git /usr/src/lancet && \
    mkdir -p /usr/src/lancet/build && cd /usr/src/lancet/build && \
    cmake -GNinja -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_BUILD_TYPE=Release .. && ninja -v && cp lancet /usr/bin/

ENTRYPOINT ["lancet"]

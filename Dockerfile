FROM gcc:10

COPY . /usr/src/lancet

RUN wget -cq "https://github.com/Kitware/CMake/releases/download/v3.18.1/cmake-3.18.1-Linux-x86_64.sh" && \
    wget -cq "https://github.com/ninja-build/ninja/releases/download/v1.10.0/ninja-linux.zip" && \
    /bin/bash "cmake-3.18.1-Linux-x86_64.sh" --skip-license --exclude-subdir --prefix=/usr && \
    apt-get update && apt-get install --yes --no-install-recommends unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && unzip -d /usr/bin ninja-linux.zip
#    git clone https://github.com/omicsnut/v2_lancet.git /usr/src/lancet && \
#    mkdir -p /usr/src/v2_lancet/build

WORKDIR /usr/src/lancet/build

RUN cmake -GNinja -DCMAKE_BUILD_TYPE=Release .. && ninja -v && cp lancet /bin/

ENTRYPOINT ["lancet"]

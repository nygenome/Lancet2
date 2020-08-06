FROM gcr.io/oss-fuzz-base/base-clang

MAINTAINER Rajeeva Musunuri <rmusunuri@nygenome.org>

RUN apt-get update && apt-get -y --no-install-recommends upgrade && \
    apt-get install -y --no-install-recommends zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev git cmake make && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p /usr/src && git clone https://github.com/omicsnut/v2_lancet.git /usr/src/lancet && \
    mkdir -p /usr/src/lancet/build && cd /usr/src/lancet/build && \
    cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Debug -DLANCET_SANITIZER=Address .. && \
    make && cp lancet /usr/bin/

ENV CC  clang
ENV CXX clang++

ENTRYPOINT ["lancet"]

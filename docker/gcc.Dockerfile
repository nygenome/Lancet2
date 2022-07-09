FROM alpine:edge

LABEL org.opencontainers.image.authors="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN apk update && apk upgrade && apk add --no-cache ca-certificates && rm -rf /var/cache/apk/* && \
    apk add --no-cache bash make cmake ninja git linux-headers libc-dev gcc g++ \
                       zlib-dev bzip2-dev xz-dev curl-dev openssl-dev

ENV CC  gcc
ENV CXX g++

RUN mkdir -p /usr/src && git clone https://github.com/nygenome/Lancet2.git /usr/src/lancet && \
    mkdir -p /usr/src/lancet/build && cd /usr/src/lancet/build && \
    cmake -GNinja -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release .. && \
    ninja -v && cd / && cp /usr/src/lancet/build/lancet /usr/bin/

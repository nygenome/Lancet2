FROM alpine:edge

MAINTAINER Rajeeva Musunuri <rmusunuri@nygenome.org>

RUN apk update && apk upgrade && apk add --no-cache ca-certificates && rm -rf /var/cache/apk/* && \
    apk add --no-cache linux-headers libc-dev gcc g++ bash make cmake \
                       ninja git zlib-dev bzip2-dev xz-dev curl-dev openssl-dev

RUN mkdir -p /usr/src && git clone https://github.com/omicsnut/v2_lancet.git /usr/src/lancet && \
    mkdir -p /usr/src/lancet/build && cd /usr/src/lancet/build && \
    cmake -GNinja -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release .. && \
    ninja -v && cp lancet /usr/bin/

ENV CC  gcc
ENV CXX g++

ENTRYPOINT ["lancet"]

# Build and run:
#   docker build --platform linux/amd64 -t clion/lancet2/dev-env:1.0 -f cpp-env-ubuntu.Dockerfile .

FROM ubuntu:20.04

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get -y install tzdata

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      build-essential \
      gcc \
      g++ \
      gdb \
      clang \
      make \
      git \
      ninja-build \
      cmake \
      autoconf \
      automake \
      locales-all \
      dos2unix \
      rsync \
      tar \
      python \
      python-dev \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl3-dev \
      libssl-dev \
  && apt-get clean \
  && apt-get purge \
  && rm -rf "/var/lib/apt/lists/*" "/tmp/*" "/var/tmp/*"

FROM python:3.12-slim
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get upgrade --yes --no-install-recommends && \
    python3 -m venv /pyenv && /pyenv/bin/pip install -U pip tqdm==4.66.2 pysam==0.22.0 interpret-core==0.5.1

FROM python:3.12-slim
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && apt-get upgrade --yes --no-install-recommends && \
    apt-get install  --yes --no-install-recommends time ca-certificates curl wget bcftools && \
    wget "https://raw.githubusercontent.com/nygenome/Lancet2/5e8731691e88cb41f74259aeb7fbbdbafd030cd8/python/filter_variants.py" && \
    wget "https://storage.googleapis.com/lancet-ml-models/somatic_ebm.lancet_6ef7ba445a.v1.pkl" && \
    python3 -m venv /pyenv && /pyenv/bin/pip install -U pip numpy==1.26.4 tqdm==4.66.2 pysam==0.22.0 interpret-core==0.5.1


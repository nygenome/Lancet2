# Scoring Somatic Variants

## Setup & Installation

The `score_variants.py` python script requires `Python 3.x` and some additional dependencies
which are ideally installed using a [virtual environment](https://docs.python.org/3/library/venv.html).

```bash
python3 -m venv --upgrade-deps pyenv
./pyenv/bin/pip install numpy==1.26.4 tqdm==4.66.2 pysam==0.22.0 interpret-core==0.5.1
```

The explainable somatic machine learning model ([`somatic_ebm.lancet_6ef7ba445a.v1.pkl`](https://storage.googleapis.com/lancet-ml-models/somatic_ebm.lancet_6ef7ba445a.v1.pkl))
is also needed to run the `score_variants.py` script.

## Usage

```bash
./pyenv/bin/python3 score_variants.py \
    lancet2_output.vcf.gz somatic_ebm.lancet_6ef7ba445a.v1.pkl \
    > lancet2_output.somatic_scoring.vcf
```

The `PASS` somatic variants can then be filtered from the scored VCF as follows.

```bash
bcftools view -f PASS -Oz -o lancet2_output.somatic_scoring.PASS.vcf.gz \
    lancet2_output.somatic_scoring.vcf
```

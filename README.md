# Background selection with local mutation rate variation in humans

This repository contains analyses and data for our paper on background
selection (BGS) and mutation rate variation, preprinted here: XXX

## Contents

`config` contains configuration files for the Snakemake pipeline.

`data` contains data used to compute expected diversity reductions due to BGS
and expected/observed diversity. See [data/readme.md] for more details.

`local_models` contains scripts used to compute expected BGS at a fine genomic
resolution (1kb).

`models` contains final predictions at 10kb coarseness. These are scaled-up
versions of the models in [local\_models/]. This directory also contains some
tables summarizing predictions and goodness-of-fit.

`tools` contains R scripts for downstream analysis, Python notebooks for
summarizing expectations, and the Snakemake pipeline for assembling input data.

## Requirements

### Python
To install Python dependencies, you can create a virtual environment and use
`pip`:

```
python -m venv my_env
source my_env/bin/activate
pip install numpy scipy pandas matplotlib jupyter-notebook
```

The `bgshr` package is used to predict the intensity of BGS and can
be installed from Github with:

```
pip install git+https://github.com/apragsdale/bgshr
```

## Reproducing data and predictions





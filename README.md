# Background selection with local mutation rate variation in humans

This repository contains analyses and data for our paper on background
selection (BGS) with mutation rate variation, preprinted [here](https://www.biorxiv.org/content/10.64898/2026.06.02.727906v2.abstract)

## Contents

[config/](config) contains configuration files for the Snakemake pipeline.

[data/](data) contains data used to compute expected diversity reductions due to BGS
and expected/observed diversity. See the `README.txt` files in its subdirectories
for further details.

[local\_models/](local_models) contains scripts used to compute expected BGS at
a fine genomic resolution (1kb).

[models/](models) contains final predictions at the 10kb scale. These are
scaled-up versions of the models in [local\_models/].
Fine-scale `B`-maps for a subset of models can be found in 
[models/B\_maps\_1kb/](models/B\_maps_1kb).
This directory also contains some tables summarizing predictions and 
goodness-of-fit.

[tools/](tools) contains R scripts for downstream analysis, Python notebooks
for summarizing expectations, and the Snakemake pipeline for assembling input 
data. See its [README.md](tools/README.md) for more.

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

To reproduce the empirical dataset, follow the instructions for running the
[Snakefile](tools/Snakefile) in [tools/README.md](tools/README.md).
Some large files must be downloaded manually, and the configuration file
[config/config.yaml](config/config.yaml) edited to point to the directory that
contains them.

Scripts for fitting the effective population size and predicting `B` are
written to run on the UW-Madison Center for High-Throughput Computing (CHTC)
server, which uses the HTCondor environment.


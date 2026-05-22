
## Preparing input data

Input data can be reproduced using the `Snakefile` in this directory.

## Building tables for downstream analysis

Because predicting B is computationally expensive, we ran predictions on a
high-throughput computing cluster.
There is therefore a discontinuity in the pipeline.
Scripts for reproducing B predictions can be found in `../local_models/`.
After predicted tables have been assembled in that directory, Snakemake can be
invoked to build tables for downstream analysis using

```
snakemake --cores 8 build_models
```

replacing `8` with the number of cores you would like to use.
This will append an `avg_pi` column containing observed pi to the predicted
tables and scale them to 10kb, so that they use less disk space.

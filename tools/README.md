
## Contents

## Running pipelines

### Input data not included in the Snakefile

Functional/phylogenetic conservation annotations within
[../data/annotations](../data/annotations) are not included in the Snakefile,
but are hosted in this repository.

### Preparing input data

Input data can be reproduced using routines of the `Snakefile` in this
directory.
These are modular, but should be run in the order below due to dependencies
between the output files.

#### Mutation rate data
To obtain mutation data, download files chromosome-specific maps 
`*_rate_v5.2_TFBS_correction_all.vcf.gz` from
https://genetics.bwh.harvard.edu/downloads/Vova/Roulette/ (~71GB).
These are VCF files which contain predicted mutation rates from the Roulette
(Seplyarskiy et al., 2023), Carlson (Carlson et al., 2018), and gnomAD
(Karczewski et al., 2020) mutation models.
In this repository, we stylize the names of the respective models as
`roulette`, `carlson` and `gnomad`.
You can place these in [../data/roulette\_vcfs/](../data/roulette_vcfs) or
elsewhere.
If you place them in a different directory, you should edit the [
config.yaml](../config/config.yaml) file to point `roulette_vcf_dir` to that
directory.
To convert mutation maps into site-resolution binary arrays (`.npy`), run:

```
snakemake --cores 4 parse_mutation_maps
```

You may wish to run this command with relatively few cores, as the process
is memory-intensive.

#### Genetic masks

We use genetic masks to determine which sites should be included in
observations, predictions, and inference.
These incorporate the strict 1000 Genomes mask, which specifies regions with
high sequence uniqueness and mapping quality, and the coverage of each
mutation map.
The `roulette` and `gnomad` maps have the same coverage, so only a single mask
labelled `roulette` is produced for these two models.
You can construct the genetic masks using:

```
snakemake --cores 8 build_genetic_masks
```

This parses the coverage of mutation maps with unique coverage profiles and
intersects covered intervals with 1KG strict mask intervals.

#### Diversity data

We estimate diversity from sample-specific GVCF files.
To reproduce the diversity dataset, download the files corresponding to 108
YRI samples (listed in [YRI\_sample\_ids.txt](../config/YRI_sample_ids.txt)) from
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data\_collections/1000G\_2504\_high\_coverage/working/20190425\_NYGC\_GATK/raw\_calls\_updated/
and update the `gvcf_dir` field in [config.yaml](../config/config.yaml).
You'll also need to update the `allele_count_dir` field in the configuration
file.
It may be desirable to place the GVCF and allele count array directories on
an external drive, as these take a large amount of disk space (~500GB each).
After gathering these files, you can regenerate diversity data by running:

```
snakemake --cores 4 estimate_diversity
```

### Computing expected maps

Because predicting `B` is computationally expensive, we ran predictions on a
high-throughput computing cluster.
There is therefore a discontinuity in the Snakemake pipeline.
Tables with predicted `B`-values can be reproduced using the HTCondor scripts
housed in [local\_models/](../local_models).

### Building tables for downstream analysis

After predicted tables have been assembled in that directory, Snakemake can be
invoked to build tables for downstream analysis using

```
snakemake --cores 8 get_model_tables
```

replacing `8` with the number of cores you would like to use.
This will append an `avg_pi` column containing observed pi to the predicted
tables and scale them to 10kb, so that they use less disk space.

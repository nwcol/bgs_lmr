## Scripts

This directory holds Python scripts invoked by the Snakemake file in the parent
directory. These are:

### Calculating observed diversity

`parse_gvcf.py`: Parses tables of reference/alternate allele counts from sample-
specific gVCF files (genome-wide VCF; these contain invariant as well as
polymorphic sites) and saves them as .npy files.

`sum_allele_count_arrays.py`: Sums sample-specific allele count arrays together
and saves the output in another .npy file.

`build_ndns_array.py`: Computes the numbers of same/different-by-state allele
copy pairs (nD/nS) at each site using multi-sample tables of reference/
alternate allele counts and saves the result in a .npy file. This is the
variation data to which we fit Ne and which we use to estimate pairwise
nucleotide diversity.

`build_pi_table.py`: Calculates diversity in regular genomic windows using an 
nD/nS table and saves the result as a .csv file.

### Processing mutation site-resolution mutation maps and map coverages

`parse_roulette_vcf.py`: Parses site-resolution mutation rates from a .vcf file
and saves the resulting array in a .npy file, with `nan` where mutation rate
data is missing.

`parse_mutation_map_coverage.py`: Creates a BED file recording intervals where
mutation rate data is non-missing, from a mutation rate .npy file.

`merge_bedfiles.py`: Merges two or more BED mask files, by either taking the
union of their intervals or the intersection. When the intersection is taken,
only sites which occur in intervals in every input BED file are represented in
intervals of the output file.

### Processing tables of B and pi predictions

`build_B_table.py`: Appends an observed diversity column to a table of B and
diversity predictions, drops windows in the centomere (and HLA locus, for
chromosome 6), and scales the table to 10kb window size. To perform scaling,
we weight observations/predictions by the number of accessible sites in each
input window (these are typically 1kb), `num_sites`. 

`build_B_map.py`: Creates a minimal table holding the chromosome number, focal
site positions, and predicted B-values:

```
chrom,pos,B
```

See the scripts themselves for more detailed documentation.

### Table format

Throughout, files referred to as `tables` are .csv files with the BEDGRAPH
format, e.g., structured as

```
chrom,chromStart,chromEnd,data_0,data_1,...
1,0,1000,0.0,0.0,...
```

To scale tables to coarser window sizes, we take averages weighted by the
number of accessible sites in each window. Therefore tables should minimally
have a `num_sites` column;

```
chrom.chromStart,chromEnd,num_sites,...
1,0,1000,2,...
1,1000,2000,450,...
```

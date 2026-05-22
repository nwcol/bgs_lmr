#!/bin/bash

# special arguments for retrieving large files
chrom=$1
mutmodel=$2
cp /staging/valadaresbar/mutation_arrays/${mutmodel}/${mutmodel}_map_chr${chrom}.npy  .
cp /staging/valadaresbar/ndns_arrays/ndns_YRI_chr${chrom}.npy .

bgshr fit_Ne ${@:3}

rm ${mutmodel}_map_chr${chrom}.npy ndns_YRI_chr${chrom}.npy
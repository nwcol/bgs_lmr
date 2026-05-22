#!/bin/bash

# special arguments for retrieving large files
chrom=$1
mutmodel=$2
cp /staging/nwcollier/mutation_arrays/${mutmodel}/${mutmodel}_map_chr${chrom}.npy .

bgshr predict_B ${@:3}

rm ${mutmodel}_map_chr${chrom}.npy

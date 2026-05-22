#!/bin/bash

chrom=$1
mutmodel=$2
cp /staging/nwcollier/mutation_arrays/${mutmodel}/${mutfile} .
cp /staging/nwcollier/ndns_arrays/ndns_YRI_chr${chrom}.npy .

bgshr fit_Ne ${@:3}

rm ${mutfile} ndns_YRI_chr${chrom}.npy
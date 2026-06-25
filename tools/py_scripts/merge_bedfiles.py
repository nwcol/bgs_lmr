"""
Merges two or more BED files. If executed with the --union flag, merges by
taking the union of intervals in input files; if executed with the --intersect
flag, merges by taking the intersection of intervals- only sites which are
covered by every input BED file will appear in the output file.

Inputs/outputs may be unzipped or gzipped.

Usage
-----
$ python merge_bedfiles.py \
   --union \
   -i input1.bed input2.bed.gz \
    -o merged_file.bed.gz

Notes
-----
- Will raise an error if neither --union nor --intersect is given.

- Chromosome numbers are checked for consistency. `chr22` and `22` are
  consistent.

- Chromosome number formats in output files will match the chromosome number 
  format of the first input file (if the first input file has `chr22`, so will
  the output).
"""

import argparse
import bgshr
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fnames",
        type=str,
        nargs="*",
        required=True,
        help="input BED files (.bed or .bed.gz)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="input BED files (.bed or .bed.gz)")
    parser.add_argument(
        "--union",
        action="store_true",
        help="take the union of intervals in input BED files")
    parser.add_argument(
        "--intersect",
        action="store_true",
        help="take the intersection of intervals in input BED files")
    return parser.parse_args()


def merge_bedfiles(in_fnames, out_fname, union=False, intersect=False):
    if (not union and not intersect) or (union and intersect):
        raise ValueError("you must use either --union or --intersect")
    regions_list = []
    chroms = []
    for in_fname in in_fnames:
        regions, chrom = bgshr.Util.read_bedfile(in_fname, get_chrom=True)
        regions_list.append(regions)
        chroms.append(chrom)
    assert len(set([c.lstrip("chr") for c in chroms])) == 1
    chrom = chroms[0]
    if intersect:
        out_regions = bgshr.Util.intersect_elements(regions_list)
    else:
        out_regions = bgshr.Util.merge_elements(regions_list)
    bgshr.Util.write_bedfile(out_fname, out_regions, chrom)
    return


def main():
    args = get_args()
    merge_bedfiles(
        args.in_fnames,
        args.out_fname,
        union=args.union,
        intersect=args.intersect)
    return


if __name__ == "__main__":
    main()


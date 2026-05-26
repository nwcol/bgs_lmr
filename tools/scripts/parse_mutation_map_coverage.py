"""
Load a .npy file representing a site-resolution mutation map, and construct a
.bed mask file recording intervals of non-missing mutation data.

Usage
-----
$ python parse_mutation_map_coverage.py \
    -i mut_map.npy \
    -chrom 22 \
    -o mut_map_coverage.bed.gz

Notes
-----
- `nan` entries in the input array are taken to be missing, and all other
  elements as non-missing.

- The `--chrom` argument exists to provide a chromosome name for the output
  BED file, as this information is not encoded in the input array.
"""

import argparse
import numpy as np
import bgshr


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fname",
        type=str,
        required=True,
        help="Pathname of input site-resolution mutation map (.npy)")
    parser.add_argument(
        "-chrom",
        "--chrom",
        type=str,
        default="0",
        help="String with which to fill BED 'chrom' column (default '0')")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="Pathname of output mask file (.bed)")
    return parser.parse_args()


def parse_mutation_map_coverage(in_fname, out_fname, chrom):
    site_map = np.load(in_fname)
    bool_mask = np.isnan(site_map)
    intervals = bgshr.Util.mask_to_elements(bool_mask)
    bgshr.Util.write_bedfile(out_fname, intervals, chrom)
    return


def main():
    args = get_args()
    parse_mutation_map_coverage(args.in_fname, args.out_fname, args.chrom)
    return


if __name__ == "__main__":
    main()


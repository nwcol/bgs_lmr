"""
Sums over many sample-specific allele count arrays. These are arrays with 
shape (L, 2), where L is chromosome length, which old counts of reference/
alternate allele copies.

Note that the data type used for these arrays (np.uint8) cannot support allele
counts > 2 ** 8; in this work, our sample consists of 108 diploids or at most
216 allele copies.

Usage:
python sum_allele_count_arrays.py -i sample0.npy sample1.npy ... \
    -o all_allele_counts.npy
"""

import argparse
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fnames",
        type=str,
        nargs="*",
        required=True,
        help="paths to input array files (.npy)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="path to output file")
    return parser.parse_args()


def sum_allele_count_arrays(in_fnames, out_fname):
    sum_arrs = 0
    for in_fname in in_fnames:
        arr = np.load(in_fname)
        summed_counts += sample_counts
    np.save(out_fname, summed_counts)
    return


def main():
    args = get_args()
    sum_allele_count_arrays(args.in_fnames, args.out_fname)
    return


if __name__ == "__main__":
    main()


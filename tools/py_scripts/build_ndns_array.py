"""
Compute the number of same/different-by-state allele copy pairs at each site,
given an array of reference/alternate allele counts at each site.

Usage
-----
$ python build_ndns_array.py -i allele_counts.npy -o ndns_array.npy

Notes
-----
- Where a site has 0 observed alleles, the output nd/ns array has 0 in each
  column.
"""

import argparse
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fname",
        type=str,
        required=True,
        help="input allele counts array (.npy)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="output filename (.npy)")
    return parser.parse_args()


def compute_ndns(nrna):
    """
    Computes an nDnS array from nRnA, the numbers of observed
    reference and alternate allele copies.

    :param nrna: Array of reference/alternate allele counts.
    """
    nrna = np.array(nrna, np.uint16)
    nr = nrna[:, 0]
    na = nrna[:, 1]
    # Number of distinct heterozygous pairs
    nd = nr * na
    n = nr + na
    nt = n * (n - 1) // 2
    ns = nt - nd
    ndns = np.stack([nd, ns], axis=1)
    return ndns


def build_ndns_array(in_fname, out_fname):
    nrna = np.load(in_fname)
    ndns = compute_ndns(nrna)
    np.save(out_fname, ndns)
    return


def main():
    args = get_args()
    build_ndns_array(args.in_fname, args.out_fname)
    return


if __name__ == "__main__":
    main()


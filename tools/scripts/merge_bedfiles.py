"""
Merges two or more BED files. If executed with the --union flag, merges by
taking the union of intervals in input files; if executed with the --intersect
flag, merges by taking the intersection of intervals- only sites which are
covered by every input BED file will appear in the output file.

Inputs/outputs may be unzipped or gzipped.

Example:
$ python merge_bedfiles.py -i input1.bed input2.bed.gz \
    -o merged_file.bed.gz --union

Will raise an error if neither --union nor --intersect is given.
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
    """
    """
    if not union and not intersect:
        raise ValueError("you must use either --union or --intersect")
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

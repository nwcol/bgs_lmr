"""
Extracts the `B` column from a table of predictions and saves it alongside the
corresponding focal sites `pos` in a .csv file.

`pos` are taken to be the centers of the windows in the input table. This is
the default setup for predictions made with ``bgshr``.

Usage:
$ python build_B_map.py -i model_predictions.csv.gz -o B_map.csv.gz
"""

import argparse
import bgshr
import numpy as np
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fname",
        type=str,
        required=True,
        help="path to input prediction table (.csv or .csv.gz)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="path for output file (.csv or .csv.gz)")
    return parser.parse_args()


def build_B_map(in_fname, out_fname):
    """
    Loads a prediction table and saves the predicted B-values/focal sites in
    a compact output table.

    :param in_fname: Path to input prediction table.
    :param out_fname: Path at which to save output.
    """
    df = pandas.read_csv(in_fname)
    # Recover focal sites where B was predicted
    starts = np.array(df["chromStart"])
    ends = np.array(df["chromEnd"])
    pos = ((starts + ends) / 2).astype(np.int64)
    df_out = pandas.DataFrame(
        {"chrom": df["chrom"], "pos": pos, "B": df["B"]})
    if out_fname.endswith(".gz"):
        df_out.to_csv(out_fname, index=False, compression="gzip")
    else:
        df_out.to_csv(out_fname, index=False)
    return


def main():
    args = get_args()
    build_B_map(args.in_fname, args.out_fname)
    return


if __name__ == "__main__":
    main()


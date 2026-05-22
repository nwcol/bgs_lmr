"""
Combines a table holding predicted B, pi with a table containing observed pi.
The goal is to create a comprehensive table of expectations/observations with
an efficient window size for downstream analysis.

Also scales the window size to 10kb; drops any empty windows (windows with 
`num_sites` == 0) after scaling; and drops windows within the centromere, as
annotated in a required chromosome band file. If the chromosome number is 6,
further drops windows within the annotated HLA region.

Usage:
$ python build_B_table.py -i model_predictions.csv.gz -pi pi_tbl.csv.gz \
    -b cytoBand.txt -o B_tbl.csv.gz
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
        "-pi",
        "--pi_fname",
        type=str,
        required=True,
        help="path to input observed pi table (.csv or .csv.gz)")
    parser.add_argument(
        "-b",
        "--band_fname",
        type=str,
        required=True,
        help="path to input chromosome band file (tab-separated text file)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="path for output file (.csv or .csv.gz)")
    return parser.parse_args()


def get_centromere(band_fname, chrom):
    """
    Get upper and lower bounds of the `acen` annotation from
    cytoBand.txt.
    """
    band_df = pandas.read_csv("../data/cytoBand.txt", sep="\t",
        names=["chrom", "start", "end", "band", "stain"])
    sub_df = band_df[band_df["chrom"] == f"chr{chrom}"]
    start = np.min(sub_df[sub_df["stain"] == "acen"]["start"])
    end = np.max(sub_df[sub_df["stain"] == "acen"]["end"])
    return (start, end)


def build_B_table(in_fname, pi_fname, band_fname, out_fname, scale=1e4):
    """
    Adds the `avg_pi` column loaded from `pi_fname` to a table of predictions
    (`in_fname`), scales to `scale`, and saves the result.

    :param in_fname: Path to input prediction table.
    :param pi_fname: Path to input diversity table.
    :param band_fname: Path to chromosome band file. Expected to lack a 
        header and contain these columns;
            chrom start end band stain
        regions where `stain` is `acen` are dropped.
    :param out_fname: Path for output file.
    :param scale: Window size to which to scale, in bp (default 1e4; 10kb).
    """
    # Load inputs
    df = pandas.read_csv(in_fname)
    df_pi = pandas.read_csv(pi_fname)

    # Scale inputs to desired `scale`
    df = bgshr.Util.scale_genome_table(df, scale)
    df_pi = bgshr.Util.scale_genome_table(df_pi, scale)

    # Drop any trailing rows from the diversity table
    if len(df_pi) > len(df):
        df_pi = df_pi.loc[:len(df) - 1]

    assert len(df_pi["avg_pi"]) == len(df["exp_pi"])
    assert np.all(np.array(df_pi["chromStart"]) == np.array(df["chromStart"]))
    assert np.all(np.array(df_pi["num_sites"]) == np.array(df["num_sites"]))

    df["avg_pi"] = df_pi["avg_pi"]

    # Drop windows within the centromere
    chrom = next(iter(df["chrom"]))
    cen_start, cen_end = get_centromere(band_fname, chrom)
    df = df[np.logical_not((df["chromStart"] >= cen_start)
                            & (df["chromEnd"] <= cen_end))]

    # Drop the HLA region from chromosome 6
    if chrom == 6:
        # determined with UCSC genome browser
        hla_start, hla_end = 28500000, 33500000
        df = df[np.logical_not((df["chromStart"] >= hla_start)
                                & (df["chromEnd"] <= hla_end))]

    # Re-order columns and drop any unneeded ones
    df = df.loc[:, [
        "chrom",
        "chromStart",
        "chromEnd",
        "num_sites",
        "del_sites",
        "avg_mut",
        "del_mut",
        "avg_rec",
        "B",
        "exp_del_pi",
        "exp_pi",
        "avg_pi"
    ]]
    # Drop empty windows
    df = df[df["num_sites"] > 0]

    if out_fname.endswith(".gz"):
        df.to_csv(out_fname, index=False, na_rep=0.0, compression="gzip")
    else:
        df.to_csv(out_fname, index=False, na_rep=0.0)
    return


def main():
    args = get_args()
    build_B_table(
        args.in_fname,
        args.pi_fname,
        args.band_fname,
        args.out_fname)
    return


if __name__ == "__main__":
    main()


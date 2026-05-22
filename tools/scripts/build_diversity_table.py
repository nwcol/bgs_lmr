"""
Computes diversity in genomic windows from an nd/ns table. Output will be 
gzipped if the output path given ends in `.gz`.

Usage:
python build_diversity_table.py -i ndns_tbl.npy --scale 1e3 \
    --mask mask_file.bed.gz -o pi_tbl.csv.gz
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
        help="path to input ndns array file (.npy)")
    parser.add_argument(
        "--scale",
        type=int,
        default=1000,
        help="window size for pi table, in bp (default 1000)")
    parser.add_argument(
        "--mask",
        type=str,
        required=True,
        help="BED file to filter sites (.bed or .bed.gz)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str,
        required=True,
        help="path for output file (.csv or .csv.gz)")
    return parser.parse_args()


def build_diversity_table(in_fname, mask_fname, out_fname, scale):
    """
    Computes diversity in regular windows, following the application of a
    genetic mask, and saves it in a .csv file.

    :param in_fname: Path to input nD,nS file (.npy)
    :param mask_fname: Path to input mask file (.bed or .bed.gz)
    :param out_fname: Path for output file (.csv or .csv.gz)
    :param scale: Window size to use, in base pairs (default 1000)
    """
    ndns = np.load(in_fname)
    L = len(ndns)
    mask_regions, chrom = bgshr.Util.read_bedfile(mask_fname, get_chrom=True)
    mask = bgshr.Util.elements_to_mask(mask_regions, L=L)
    # Diversity equals nD / (nD + nS)
    nD = ndns[:, 0]
    nS = ndns[:, 1]
    nT = nD + nS
    covered = nT > 0
    site_pi = np.zeros(L)
    site_pi[covered] = nD[covered] / nT[covered]
    site_pi = np.ma.array(site_pi, mask=mask)
    windows = np.stack((np.arange(0, L - scale, scale, dtype=np.int64),
                        np.arange(scale, L, scale, dtype=np.int64)), axis=1)
    avg_pi, num_sites = bgshr.Util.compute_window_averages(windows, site_pi)
    data = {
        "chrom": [chrom] * len(windows),
        "chromStart": windows[:, 0],
        "chromEnd": windows[:, 1],
        "num_sites": num_sites, 
        "avg_pi": avg_pi}
    df = pandas.DataFrame(data)
    if out_fname.endswith(".gz"):
        df.to_csv(out_fname, index=False, compression="gzip")
    else:
        df.to_csv(out_fname, index=False)
    return


def main():
    args = get_args()
    build_diversity_table(
        args.in_fname,
        args.mask,
        args.out_fname,
        args.scale)
    return


if __name__ == "__main__":
    main()


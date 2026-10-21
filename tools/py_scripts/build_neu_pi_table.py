"""
Compute neutral diversity in regular genomic windows from an nd/ns table.
Neutral sites are sites which:
    1. Pass an accessibility mask (--mask)
    2. Are not in *any* of a set of masks specifying constrained sequence
       (--elements)
"""

import argparse
import bgshr
import numpy as np
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
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
        "--elements",
        type=str,
        nargs="*",
        help="BED file(s) specifying sites under selective constraint"
    )
    parser.add_argument(
        "--drop",
        type=str,
        nargs="*",
        help="If given, drop windows with 0 total sites"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="path for output file (.csv or .csv.gz)")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    ndns = np.load(args.input)
    L = len(ndns)
    # Load mask
    base_regions, chrom = bgshr.Util.read_bedfile(args.mask, get_chrom=True)
    # Load selectively constrained elements
    elements = [bgshr.Util.read_bedfile(f) for f in args.elements]
    # No need to resolve elements (purge overlaps) because we are going to add
    # them (take the union of their sites) anyway.
    comb_elements = bgshr.Util.merge_elements(elements)
    # Excise constrained sites from the mask
    mask_regions = bgshr.Util.subtract_elements(base_regions, comb_elements)
    mask = bgshr.Util.elements_to_mask(mask_regions, L=L)
    # Diversity equals nD / (nD + nS)
    nD = ndns[:, 0]
    nS = ndns[:, 1]
    nT = nD + nS
    covered = nT > 0
    site_pi = np.zeros(L)
    site_pi[covered] = nD[covered] / nT[covered]
    site_pi = np.ma.array(site_pi, mask=mask)
    windows = np.stack(
        (np.arange(0, L - args.scale, args.scale, dtype=np.int64),
         np.arange(args.scale, L, args.scale, dtype=np.int64)), axis=1)
    neu_pi, neu_sites = bgshr.Util.compute_window_averages(windows, site_pi)
    # Drop windows with 0 total sites if requested
    if args.drop:
        base_mask = bgshr.Util.elements_to_mask(base_regions, L=L)
        dummy = np.ma.ones(L, mask=base_mask)
        _, num_sites = bgshr.Util.compute_window_averages(windows, dummy)
        windows = windows[num_sites > 0]
        neu_sutes = neu_sites[num_sites > 0]
        neu_pi = neu_pi[num_sites > 0]
    data = {
        "chrom": [chrom] * len(windows),
        "chromStart": windows[:, 0],
        "chromEnd": windows[:, 1],
        "neu_sites": neu_sites,
        "neu_pi": neu_pi}
    df = pandas.DataFrame(data)
    if args.output.endswith(".gz"):
        df.to_csv(args.output, index=False, compression="gzip")
    else:
        df.to_csv(args.output, index=False)


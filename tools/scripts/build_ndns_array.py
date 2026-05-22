"""

"""

import numpy as np


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
    nrna = np.load(f"../nrna_arrays/YRI/nrna_YRI_chr{chrom}.npy")
    ndns = compute_ndns(nrna)
    np.save(f"YRI/ndns_YRI_chr{chrom}.npy", ndns)
    return


def main():

    return


if __name__ == "__main__":
    main()

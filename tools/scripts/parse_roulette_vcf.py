"""
Parse a Roulette VCF file holding mutation rates and save the vector of rates 
in a binary .npy file. The output is 0-indexed, with np.nan at sites where
mutation data is missing from the input file. The `field` argument determines
which mutation map will be parsed. Possible choices are:

field   mutation model  publication
MR      Roulette        Seplyarskiy et al 2023
MC      Carlson         Carlson et al 2018
MG      gnomAD          Karczewski et al 2020

Usage:
$ python parse_roulette_vcf.py -i input.vcf.gz -o output.npy --field MR
"""

import argparse
import gzip
import numpy as np
import re
from bgshr import Util


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--in_fname",
        type=str,
        required=True,
        help="path to input Roulette file (.vcf.gz)")
    parser.add_argument(
        "-o",
        "--out_fname",
        type=str, 
        required=True,
        help="path for output file (.npy)")
    parser.add_argument(
        "-field",
        "--field",
        type=str,
        default="MR",
        help="mutation model to parse (MR, MC or MG; see script docstring)")
    parser.add_argument(
        "-v",
        "--verbosity",
        type=int,
        default=0,
        help="prints progress messages if > 0 (default 0)")
    return parser.parse_args()


def read_sequence_length(fname):
    """
    Obtains the sequence length from the input file header.
    """
    seq_lens = dict()
    with gzip.open(fname, "rb") as file:
        for lineb in file:
            line = lineb.decode()
            if line.startswith("#"):
                if line.startswith("##contig=<ID"):
                    _, _chrom, _length, *__ = re.split("<|>|,", line)
                    chrom_num = _chrom.split("=")[1]
                    length = int(_length.split("=")[1])
                    # we only want autosomes;
                    if chrom_num.strip("chr").isnumeric():
                        seq_lens[chrom_num] = length
                else:
                    pass
            else:
                chrom_num = line.split()[0]
                break
    sequence_length = seq_lens[chrom_num]
    return sequence_length


def read_vcf(in_fname, sequence_length, field, verbosity=1e6):
    """
    Reads a vector of mutation rates from a VCF file.

    Leaves missing data as np.nan.

    :param in_fname: Path to input file. Expected to be gzipped (.gz)
    :param sequence_length: Length of the chromosome- parsed from input file
        header.
    :param field: Mutation rate annotation to access (see script docstring)
    """
    rates = np.full(sequence_length, np.nan, dtype=np.float64)
    i = 0
    with gzip.open(in_fname, "rb") as file:
        for lineb in file:
            line = lineb.decode()
            if line.startswith("#"):
                continue
            split_line = line.split()
            pos1 = int(split_line[1])
            pos0 = pos1 - 1
            _info = split_line[7]
            info = {x: y for (x, y) in [x.split("=") for x in _info.split(";")]}
            if field not in info:
                continue
            if np.isnan(rates[pos0]):
                rates[pos0] = 0
            rates[pos0] += float(info[field])
            i += 1
            if verbosity > 0 and i % verbosity == 0:
                print(Util._get_time(), f"parsed {i} lines")
    print(Util._get_time(), f"finished; parsed {i} lines")
    return rates


def parse_roulette_vcf(in_fname, out_fname, field, verbosity):
    """
    Loads mutation rates, scales them by the appropriate factor, and saves them 
    in an .npy file.

    Coefficients are from:
    https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate

    :param in_fname: Path to input file (.vcf.gz)
    :param out_fname: Path for output file (.npy)
    :param field: Mutation annotation to access (see script docstring)
    """
    if field == "MR":
        c = 1.015e-7 / 2
    elif field == "MC":
        c = 2.086e-9 / 2
    elif field == "MG":
        c = 1.015e-7 / 2
    else:
        raise ValueError("Unrecognized `--field`")
    sequence_length = read_sequence_length(in_fname)
    arr = read_vcf(in_fname, sequence_length, field, verbosity=verbosity)
    arr *= c
    np.save(out_fname, arr)
    return


def main():
    args = get_args()
    parse_roulette_vcf(
        args.in_fname,
        args.out_fname,
        args.field,
        args.verbosity)
    return


if __name__ == "__main__":
    main()


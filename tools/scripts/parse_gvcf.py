"""
Parse a table of reference/alternate allele counts for SNPs in a one-sample
gVCF file (a .vcf-format file covering both variant and invariant sites).

Outputs a table of allele counts in .npy format for each chromosome. The 0th
column of this table holds reference allele counts and the 1st column holds
alternate allele counts (alternate alleles are not distinguished).

Usage
-----
$ python parse_gvcf.py \
    -i input.vcf.gz \
    -g hg38.genome \
    -o allele_counts.npy
"""

import argparse
from dataclasses import dataclass
import gzip
import numpy as np
import re
from bgshr import Util


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--in_fname',
        type=str,
        required=True,
        help='Input gVCF file with .vcf or .vcf.gz extension'
    )
    parser.add_argument(
        '-o',
        '--out_prefix',
        type=str,
        required=True,
        help='Prefix for output filenames; writes to {out_prefix}_chrXX.npy'
    )
    parser.add_argument(
        '-g',
        '--genome_file',
        type=str,
        default=None,
        help='Genome file holding chromosome 1-22 lengths'
    )
    parser.add_argument(
        '--min_GQ',
        type=float,
        default=30,
        help='Minimum GQ to pass filter'
    )
    parser.add_argument(
        '--min_QUAL',
        type=float,
        default=50,
        help='Minimum QUAL to pass filter'
    )
    parser.add_argument(
        '--verbose',
        type=int,
        default=1000000,
        help='Line interval for printing status reports'
    )
    return parser.parse_args()


def read_seq_lens(fname):
    """
    Obtain sequence lengths from the input file header.
    """
    open_func = gzip.open if fname.endswith('gz') else open
    seq_lens = dict()
    with open_func(fname, 'rb') as file:
        for lineb in file:
            line = lineb.decode()
            if line.startswith('#'):
                if line.startswith('##contig=<ID'):
                    _, _chrom, _length, __ = re.split('<|>|,', line)
                    chrom_num = _chrom.split('=')[1]
                    length = int(_length.split('=')[1])
                    # we only want autosomes;
                    if chrom_num.strip('chr').isnumeric():
                        seq_lens[chrom_num] = length
                else:
                    pass
            else:
                break
    return seq_lens


def read_genome_file(fname):
    """
    Read autosomal chromosome lengths from a .genome file.
    """
    seq_lens = {}
    with open(fname, 'r') as file:
        for line in file:
            chrom, length = line[:2].split()
            if length.isnumeric():
                seq_lens[chrom] = length
    return seq_lens


def get_genotype_cache(phased=False):
    """
    Get a dictionary mapping genotypes to reference/alternate allele counts.
    """
    sep = '|' if phased else '/'
    is_ref = lambda x : 1 if x == 0 else 0
    cache = dict()
    for allele1 in range(4):
        for allele2 in range(4):
            genotype = f'{allele1}{sep}{allele2}'
            num_refs = is_ref(allele1) + is_ref(allele2)
            num_alts = 2 - num_refs
            allele_counts = [num_refs, num_alts]
            cache[genotype] = allele_counts
    return cache


def write_log(chrom_stats, prefix):
    """
    write a tab seperated log file with per-chromosome statistics
    """
    fname = f'{prefix}_log.txt'
    keys = list(chrom_stats['chr1'].__dict__.keys())
    header = 'chrom\t' + '\t'.join(keys) + '\n'
    with open(fname, 'w') as file:
        file.write(header)
        for chrom in chrom_stats:
            vals = [str(chrom_stats[chrom].__dict__[key]) for key in keys]
            line = str(chrom) + '\t' + '\t'.join(vals) + '\n'
            file.write(line)
    return


def main():
    args = get_args()
    fin = args.in_fname
    prefix = args.out_prefix
    min_GQ = args.min_GQ
    min_QUAL = args.min_QUAL
    if args.genome_file is None:
        seq_lens = read_seq_lens(fin)
    else:
        seq_lens = read_genome_file(args.genome_file)
    gt_cache = get_genotype_cache()
    # caches and statistics
    format_cache = dict()
    chrom_stats = dict()

    @dataclass
    class Stats:
        num_lines: int = 0
        num_hom_refs_failed_gq: int = 0
        num_hom_refs_passed: int = 0
        num_snps_failed_gq: int = 0
        num_snps_failed_qual: int = 0
        num_snps_passed: int = 0
        num_mnvs: int = 0
        len_mnvs: int = 0

    last_chrom = None
    open_func = gzip.open if fin.endswith('gz') else open
    print(Util._get_time(), f'parsing allele counts from {fin}')
    with open_func(fin, 'rb') as file:
        for lineb in file:
            line = lineb.decode()
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, qual, __, info, format, sample = \
                line.split()
            # handle chromosome switching
            if chrom != last_chrom:
                if last_chrom is not None:
                    fname = f'{prefix}_{last_chrom}'
                    np.save(fname, counts)
                    dic = stats.__dict__
                    prstats = '\t'.join(f'{key} = {dic[key]}' for key in dic)
                    print(
                        Util._get_time(), 
                        f'wrote {last_chrom} counts to {fname}'
                        f' {prstats}'
                    )
                    chrom_stats[last_chrom] = stats
                if chrom in seq_lens:
                    counts = np.zeros((seq_lens[chrom], 2), dtype=np.uint8)
                    stats = Stats()
                else:
                    break
            last_chrom = chrom
            stats.num_lines += 1

            split_sample = sample.split(':')
            passed_GQ = True
            # filter by GQ
            if format in format_cache:
                GQ_idx = format_cache[format]
            else:
                split_format = format.split(':')
                GQ_idx = split_format.index('GQ')
                format_cache[format] = GQ_idx
            if float(split_sample[GQ_idx]) < min_GQ:
                passed_GQ = False

            hom_ref = True if alt == '<NON_REF>' else False
            start = int(pos) - 1
            if hom_ref:
                end = int(info.split('=')[1])
                if passed_GQ:
                    counts[start:end, :] = [2, 0]
                    stats.num_hom_refs_passed += end - start
                else:
                    stats.num_hom_refs_failed_gq += end - start
            else:
                # check whether the site is a SNP
                len_ref = len(ref)
                alts = alt.split(',')[:-1]
                alt_lens = sum([len(v) for v in alts])
                is_snp = True if len_ref == 1 and alt_lens == len(alts) else False
                if is_snp:
                    if passed_GQ:
                        # filter by QUAL
                        if float(qual) < min_QUAL:
                            stats.num_snps_failed_qual += 1
                        else:
                            genotype = split_sample[0]
                            counts[start, :] = gt_cache[genotype]
                            stats.num_snps_passed += 1
                    else:
                        stats.num_snps_failed_gq += 1
                else:
                    stats.num_mnvs += 1
                    stats.len_mnvs += len_ref
    write_log(chrom_stats, prefix)
    return


if __name__ == '__main__':
    main()

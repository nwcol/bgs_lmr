# Variables: chrom,B,mutmodel,mutfname,outstem

import numpy as np
import pandas


# load unlinked B-values
Ne_df = pandas.read_csv("../2026_05_12_corrected_regulatory_fits/Ne_tbl.csv")
B_df = pandas.read_csv("../../data/B_unlinked_tbl.csv")

memory = [32]*16 + [16]*6

chroms = list(range(1, 23))
mut_models = ["roulette", "gnomad", "carlson"]


for cons_model in ["split_cds_regulatory", "merged_cds_regulatory"]:
    with open(f"{cons_model}_variables.txt", "w") as fout:
        for chrom in chroms:
            for mut_model in mut_models:
                mut_file = f"{mut_model}_map_chr{chrom}.npy"

                B = next(iter(B_df[
                    (B_df["chrom"] == chrom)
                    & (B_df["mut_model"] == mut_model)
                    & (B_df["cons_model"] == cons_model)
                ]["B"]))

                Ne = next(iter(Ne_df[
                    (Ne_df["chrom"] == chrom)
                    & (Ne_df["mut_model"] == mut_model)
                    & (Ne_df["cons_model"] == cons_model)
                ]["Ne"]))

                mem_req = str(memory[chrom - 1]) + "GB"
                outstem = f"{cons_model}_{mut_model}_chr{chrom}.csv"
                variables = [
                    str(chrom),
                    str(B),
                    str(Ne),
                    mut_model,
                    mem_req,
                    outstem
                ]
                fout.write(",".join(variables) + "\n")


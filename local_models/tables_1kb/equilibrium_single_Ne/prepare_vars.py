# Variables: chrom,B,Ne,mutmodel,mutfname,outstem

import numpy as np
import pandas


# load unlinked B-values and Ne
B_df = pandas.read_csv("../../models/B_unlinked_tbl.csv")
Ne_df = pandas.read_csv("../../models/equilibrium_single_Ne/avg_Ne_tbl.csv")
memory = [32]*16 + [16]*6
chroms = list(range(1, 23))
mut_models = [
    "carlson",
    "gnomad",
    "roulette"]
cons_models = [
    "merged_cds_regulatory",
    "merged_cds_phastcons",
    "split_cds_regulatory",
    "split_cds_phastcons"]

for cons_model in cons_models:
    with open(f"{cons_model}/{cons_model}_vars.txt", "w") as fout:
        for chrom in chroms:
            for mut_model in mut_models:
                mut_file = f"{mut_model}_map_chr{chrom}.npy"

                B = next(iter(B_df[
                    (B_df["chrom"] == chrom)
                    & (B_df["mut_model"] == mut_model)
                    & (B_df["cons_model"] == cons_model)
                ]["B"]))

                Ne = next(iter(Ne_df[
                    (Ne_df["mut_model"] == mut_model)
                    & (Ne_df["cons_model"] == cons_model)
                ]["avg_Ne"]))

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

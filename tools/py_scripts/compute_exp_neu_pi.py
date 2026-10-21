"""
Add a column recording predicted neutral diversity to a chromosome landscape
table that lacks this column.

Because
    E[pi_tot] = (E[pi_del] * n_del + E[pi_neu] * n_neu) / n_tot
where n_tot = n_del + n_neu, we compute pi_neu using
    E[pi_neu] = (E[pi_tot] * n_tot - E[pi_del] * n_del) / n_tot
"""

import argparse
import numpy as np
import pandas


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    df = pandas.read_csv(args.input)
    if "exp_neu_pi" in df:
        raise ValueError("input table already has `exp_neu_pi` column")
    n_tot = np.array(df["num_sites"])
    n_del = np.array(df["del_sites"])
    n_neu = n_tot - n_del
    pi_tot = np.array(df["exp_pi"])
    pi_del = np.array(df["exp_del_pi"])
    sum_pi_neu = (pi_tot * n_tot - pi_del * n_del)
    # Save divison (avoid dividing by zero)
    pi_neu = np.zeros(len(pi_tot), dtype=np.float64)
    safe = n_neu > 0
    pi_neu[safe] = sum_pi_neu[safe] / n_neu[safe]
    df["exp_neu_pi"] = pi_neu
    # Insert exp_neu_pi before exp_del_pi
    cols = list(df.columns)
    cols.pop(-1)
    to_idx = cols.index("exp_del_pi")
    cols.insert(to_idx, "exp_neu_pi")
    df = df[cols]
    df.to_csv(args.output)


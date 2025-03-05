import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from DEA import run_dea


def bootstrap_resample(df: pd.DataFrame, design: str | pd.DataFrame) -> pd.DataFrame:
    N = len(df.columns) // 2

    if isinstance(design, pd.DataFrame):
        if "Condition" not in design.columns:
            raise Exception("Custom desgin matrix must have column 'Condition'")
        dd = design.value_counts("Condition").loc[list(set(design["Condition"]))]
        if dd.shape != (2,):
            raise Exception("Must have two conditions")
        N_control = dd.iloc[0]
        N_perturbed = dd.iloc[1]

        bootstrap_samples_c = np.random.choice(df.columns[:N_control], N_control)
        bootstrap_samples_p = np.random.choice(df.columns[N_perturbed:], N_perturbed)
        bs = list(bootstrap_samples_c) + list(bootstrap_samples_p)
        df_trial = df[bs]

    elif design == "paired":
        # preserve matched samples
        ind = np.array(np.random.choice(range(0, N), N))
        ind = np.concatenate([ind, ind + N])
        ind = np.sort(ind)
        df_trial = df.iloc[:, ind]

    elif design == "unpaired":
        bootstrap_samples_c = np.random.choice(df.columns[:N], N)
        bootstrap_samples_p = np.random.choice(df.columns[N:], N)
        bs = list(bootstrap_samples_c) + list(bootstrap_samples_p)
        df_trial = df[bs]

    return df_trial


def run_trial(savepath: str, name: str, trial_number: int, count_matrix_path: str, design: str) -> None:
    np.random.seed(trial_number)

    df = pd.read_csv(count_matrix_path, index_col=0)

    created_bootstrapped_design = False
    if trial_number == 0:  # Original, unbootstrapped df
        df_trial = df
    else:
        if design in ["paired", "unpaired"]:
            if len(df.columns) % 2 != 0:
                raise Exception("Must have balanced number of replicates per condition for paired or unpaired designs")
            df_trial = bootstrap_resample(df, design)

        elif os.path.isfile(design):
            meta = pd.read_csv(design, index_col=0)
            df_trial = bootstrap_resample(df, meta)
            meta_sub = meta.loc[df_trial.columns]
            meta_sub.copy()
            design = f"{savepath}/{name}_design_trial_{trial_number}.csv"
            meta_sub.index = pd.Index([col + str(i) for i, col in enumerate(meta_sub.index)])
            meta_sub.to_csv(design)
            created_bootstrapped_design = True

        else:
            raise Exception("Invalid desing:", design)

        # Ensure no duplicate col names
        df_trial.columns = [col + str(i) for i, col in enumerate(df_trial.columns)]

    if trial_number == 0:
        outfile = Path(f"{savepath}/{name}_original.csv")
    else:
        outfile = Path(f"{savepath}/{name}_trial_{trial_number}.csv")
    run_dea(df_trial, str(outfile), "edger", True, verbose=False, lfc=0, design=design)

    if trial_number > 0:
        tab = pd.read_csv(outfile, index_col=0)
        tab["Trial"] = trial_number
        tab.to_csv(outfile)

    # Clean up
    if created_bootstrapped_design:
        os.system(f"rm {design}")


if __name__ == "__main__":
    savepath = sys.argv[1]
    name = sys.argv[2]
    trial_number = int(sys.argv[3])
    count_matrix_path = sys.argv[4]
    design = sys.argv[5]

    CREATE_DUMMY_DATA = False

    if CREATE_DUMMY_DATA:
        df = pd.DataFrame(np.random.normal(0, 1, (10, 2)), index=range(10), columns=["logFC", "FDR"])
        df["Trial"] = trial_number
        df.to_csv(f"{savepath}/{name}_trial_{trial_number}.csv")

    else:
        run_trial(savepath, name, trial_number, count_matrix_path, design)

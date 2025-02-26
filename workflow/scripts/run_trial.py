import os
import sys

import numpy as np
import pandas as pd
from pathlib import Path

from DEA import run_dea

def bootstrap_resample(df, design):

    N = len(df.columns) // 2

    # preserve matched samples
    if design == "paired":
        ind = np.array(np.random.choice(range(0, N), N))
        ind = np.concatenate([ind, ind + N])
        ind = np.sort(ind)
        df_trial = df.iloc[:, ind]

    else:
        bootstrap_samples_N = np.random.choice(df.columns[:N], N)
        bootstrap_samples_T = np.random.choice(df.columns[N:], N)
        bs = list(bootstrap_samples_N)+list(bootstrap_samples_T)
        df_trial = df[bs]

    return df_trial

def run_trial(savepath, name, trial_number, count_matrix_path, design, seed="trial"):

    if seed == "trial":
        np.random.seed(trial_number)
    elif isinstance(seed, int):
        np.random.seed(seed)

    df = pd.read_csv(count_matrix_path, index_col=0)

    # TO DO: unbalanced number of samples per condition
    if len(df.columns) % 2 != 0:
        raise Exception("Must have balanced number of replicates per condition for now")
    
    if trial_number == 0: # Original, unbootstrapped df
        df_trial = df
        created_bootstrapped_design = False
    else:
        df_trial = bootstrap_resample(df, design)
    
        if os.path.isfile(design):
            meta = pd.read_csv(design, index_col=0)
            meta_sub = meta.loc[df_trial.columns]
            meta_sub.copy()
            design = f"{savepath}/{name}_design_trial_{trial_number}.csv"
            meta_sub.index = [col+str(i) for i, col in enumerate(meta_sub.index)]
            meta_sub.to_csv(design)
            created_bootstrapped_design = True

         # Ensure no duplicate col names
        df_trial.columns = [col+str(i) for i, col in enumerate(df_trial.columns)]
    
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

    create_dummy_Data = False

    if create_dummy_Data:
        df = pd.DataFrame(np.random.normal(0,1,(10,2)),
                        index=range(10),
                        columns=["logFC","FDR"])
        df["Trial"] = trial_number
        df.to_csv(f"{savepath}/{name}_trial_{trial_number}.csv")

    else:
        run_trial(savepath, name, trial_number, count_matrix_path, design)
    
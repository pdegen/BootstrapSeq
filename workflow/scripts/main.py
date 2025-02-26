import sys

import numpy as np
import pandas as pd

def snakemake_bootstrap(savepath, name, trial_number):
    pass

if __name__ == "__main__":

    savepath = sys.argv[1]
    name = sys.argv[2]
    trial_number = sys.argv[3]

    # dummy df for now
    df = pd.DataFrame(np.random.normal(0,1,(10,2)),
                      index=range(10),
                      columns=["logFC","FDR"])
    df["Trial"] = trial_number
    df.to_csv(f"{savepath}/{name}_trial_{trial_number}.csv")

    snakemake_bootstrap(savepath, name, trial_number)
    
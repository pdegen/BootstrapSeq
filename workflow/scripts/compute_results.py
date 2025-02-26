import sys
import os
import glob
import logging
import json

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from bootstrap import compute_spearmans

def process_results(savepath, name, trials, make_figs):

    merged_trials = pd.read_csv(Path(f"{savepath}/{name}_trials_merged_{trials}.csv"), index_col=0)
    original_results = pd.read_csv(Path(f"{savepath}/{name}_original.csv"), index_col=0)
    spearmans = compute_spearmans(original_results, merged_trials)

    outfile = Path(f"{savepath}/{name}_stats.json")

    dictionary = {
        "spearman_median": np.median(spearmans),
        "spearman_mean": np.mean(spearmans),
        "spearman_std": np.std(spearmans),
        "spearmans": spearmans.tolist()
    }
    
    json_object = json.dumps(dictionary, indent=4)
    
    with open(outfile, "w") as f:
        f.write(json_object)

    if make_figs:
        plot(spearmans, savepath, name)

def plot(spearmans, savepath, name):
    # dummy plot for now
    sns.kdeplot(spearmans)
    plt.savefig(f"{savepath}/{name}_spearman.pdf")

if __name__ == "__main__":

    savepath = sys.argv[1]
    name = sys.argv[2]
    trials = sys.argv[3]
    make_figs = sys.argv[4]
  
    process_results(savepath, name, trials, make_figs)
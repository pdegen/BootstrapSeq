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

def process_results(savepath, name, trials, make_figs):

    merged_trials = pd.read_csv(Path(f"{savepath}/{name}_trials_merged_{trials}.csv"))

    outfile = Path(f"{savepath}/{name}_stats.json")

    dictionary = {
        "spearman_median": merged_trials["logFC"].median(),
        "spearman_std": merged_trials["logFC"].std(),
    }
    
    json_object = json.dumps(dictionary, indent=4)
    
    with open(outfile, "w") as f:
        f.write(json_object)

    if make_figs:
        plot(merged_trials, savepath, name)

def plot(merged_trials, savepath, name):
    sns.kdeplot(merged_trials["logFC"])
    plt.savefig(f"{savepath}/{name}_spearman.pdf")

if __name__ == "__main__":

    savepath = sys.argv[1]
    name = sys.argv[2]
    trials = sys.argv[3]
    make_figs = sys.argv[4]
  
    process_results(savepath, name, trials, make_figs)
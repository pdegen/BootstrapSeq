import os
import re
import glob
import logging
import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from DEA import run_dea

def compute_spearmans(tab_reference, bootstrap_results):
    spearmans = []
    trials = set(bootstrap_results["Trial"])
    genes = len(bootstrap_results)//len(trials)

    # DESeq2 logFC can return nan
    tab_reference = tab_reference["logFC"].dropna()

    if len(bootstrap_results) % genes != 0:
        logging.warning("Unequal lengths, spearman not computed")
        return 
    
    for trial in sorted(trials):
        boot = bootstrap_results[bootstrap_results["Trial"]==trial]["logFC"].dropna()
        common = tab_reference.index.intersection(boot.index)
        tab_reference_common = tab_reference.loc[common].rank()
        boot = boot.loc[common].rank()
        spearman = spearmanr(tab_reference_common, boot).statistic
        spearmans.append(spearman)
    return np.array(spearmans)

def open_bootstrap_results(save_path, 
                           method, 
                           name, 
                           return_df=True
                           ):
    
    results_file = f"{save_path}/{name}.boot.trials*.{method}.csv"
    matched_files = glob.glob(results_file)
    if matched_files:
        results_file = matched_files[0]
        existing_trials = int(results_file.split(".trials")[1].split(".")[0])
        if return_df and existing_trials > 0:
            return pd.read_csv(results_file, index_col=0), results_file, existing_trials
        else:
            return None, results_file, existing_trials
            
    logging.info(f"No bootstrap results file found: {save_path}")
    results_file = f"{save_path}/{name}.boot.trials0.{method}.csv"
    return None, results_file, 0

def bootstrap_data(df, 
                   save_path, 
                   lfc, 
                   design, 
                   method, 
                   name,  
                   trials, 
                   meta=None, 
                   logfile=None,
                   maxiter=5
                   ):

    results = None
    
    # TO DO: unbalanced number of samples per condition
    if len(df.columns) % 2 != 0:
        raise Exception("Must have balanced number of replicates per condition for now")
    
    N = len(df.columns) // 2
    os.system(f"mkdir -p {save_path}/tmp")

    # If exists, read results file where trial results will be concatenated to
    results, results_file, existing_trials = open_bootstrap_results(save_path, method, name)
    
    if results is None:
        logging.info("Initializing resultsfile")
        results_file = results_file.replace('trials*','trials0')
        os.system(f"touch {results_file}")
        existing_trials = 0

    if existing_trials >= trials :
        logging.info(f"Already have {existing_trials} trials, returning")
        return "returned_early"
            
    for trial in range(existing_trials+1, trials+1):

        outfile_dea = f"{save_path}/tmp/tab.tmp.trial{trial}.csv"
        
        np.random.seed(trial) # important for multiprocessing
        print(maxiter)
        # max 5 attempts if DEA fails for small N (matrix not full rank error if too many covariates)
        for a in range(1, maxiter): 
            print(a)
            # preserve matched samples
            if design == "paired":
                ind = np.array(np.random.choice(range(0, N), N))
                ind = np.concatenate([ind, ind + N])
                ind = np.sort(ind)
                df_bag = df.iloc[:, ind]

            else:
                bootstrap_samples_N = np.random.choice(df.columns[:N], N)
                bootstrap_samples_T = np.random.choice(df.columns[N:], N)
                bs = list(bootstrap_samples_N)+list(bootstrap_samples_T)
                df_bag = df[bs]

            logging.info(f"Running trial: {trial}, samples: {df_bag.columns}, path: {save_path}")
    
            if design == "custom":
                meta_sub = meta.loc[df_bag.columns]
                meta_sub.copy()
                # add suffix to filename avoid multiprocess conflict
                design_sub = f"{save_path}/tmp/design.trial{trial}.csv"
                meta_sub.index = [col+str(i) for i, col in enumerate(meta_sub.index)]
                meta_sub.to_csv(design_sub)
            elif design in ["paired", "unpaired"]:
                design_sub = design
        
            df_bag.columns = [col+str(i) for i, col in enumerate(df_bag.columns)]
    
            run_dea(df_bag, str(outfile_dea), method, True, verbose=False, lfc=lfc, design=design_sub)
            if a > 2 and logfile is not None:
                log = f"{save_path} {name} attempts: {a}"
                os.system(f"echo {log} >> {logfile}")
            break


        trial_results = pd.read_csv(outfile_dea, index_col=0)
        trial_results["Trial"] = trial

        if results is None:
            results = trial_results            
        else:
            results = pd.concat([results,trial_results])
            

        # Increment file name, save new bootstrap result, rm old
        # TO DO: append instead
        results_file_p1 = re.sub(r"trials(\d+)", lambda m: f"trials{int(m.group(1)) + 1}", results_file)
        results.to_csv(results_file_p1)
        os.system(f"rm {results_file}")
        results_file = results_file_p1

    if logfile is not None:
        now = datetime.datetime.now()
        log = f"{save_path} {name} trials: {trials} {now}"
        os.system(f"echo {log} >> {logfile}")
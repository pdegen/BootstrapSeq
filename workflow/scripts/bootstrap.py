import datetime
import glob
import logging
import os
import re
from typing import Optional
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from DEA import run_dea


def compute_spearmans(tab_reference: pd.DataFrame | pd.Series, merged_trials: pd.DataFrame) -> Optional[np.ndarray]:
    """Compute logFC Spearman rank correlation for each trial relative to a reference.

    Parameters
    ----------
    tab_reference : pandas.DataFrame
        Output table from edegR
    merged_trials : pandas.DataFrame
        Output table from bootstrap_data()

    Returns
    -------
    numpy.array
        1D array of Spearman correlations for each trial.
    """
    spearmans = []
    trials = set(merged_trials["Trial"])
    genes = len(merged_trials) // len(trials)

    # DESeq2 logFC can return nan
    tab_reference = tab_reference["logFC"].dropna()

    if len(merged_trials) % genes != 0:
        logging.warning("Unequal lengths, spearman not computed")
        return None

    for trial in sorted(trials):
        boot = merged_trials[merged_trials["Trial"] == trial]["logFC"].dropna()
        common = tab_reference.index.intersection(boot.index)
        tab_reference_common = tab_reference.loc[common].rank()
        boot = boot.loc[common].rank()
        spearman, _ = spearmanr(tab_reference_common, boot)
        if isinstance(spearman, (int, float)) and not np.isnan(spearman):
            spearmans.append(spearman)
    return np.array(spearmans)


def open_bootstrap_results(
    save_path: str, method: str, name: str, return_df: bool = True
) -> Tuple[Optional[pd.DataFrame], str, int]:
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


def bootstrap_data(
    df: pd.DataFrame,
    save_path: str,
    lfc: float,
    design: str,
    method: str,
    name: str,
    trials: int,
    meta: Optional[pd.DataFrame] = None,
    logfile: Optional[str] = None,
    maxiter: int = 1,
):
    """Repeatedly estimate logFC on bootstrapped resamples of df using edegR or DESeq2. Stores output in merged csv
    table.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe with raw read counts.
    save_path : str
        Path to save results in.
    lfc : float
        Formal log2 fold change threshold for differential expression testing.
    design : str
        Either "custom", "paired", or "unpaired". If "custom", must provide meta. Else, df must have even number of
        columns, sorted by condition and sample number.
    method : str
        Either "edegR" or "deseq2".
    name : str
        String to tag filenames with.
    trials : int
        _description_
    meta : pandas.DataFrame, optional
        Dataframe with covariates. Index must match df columns. Must have a column called "Condition". By default None
    logfile : str, optional
        Path to logfile, by default None
    maxiter : int, optional
        How many times to re-attempt differential expression analysis with new resamples in case of failure,
        by default 1

    Raises
    ------
    Exception
        Provided df has unequal number of replicates per condition.
    """
    results = None

    # TO DO: unbalanced number of samples per condition
    if len(df.columns) % 2 != 0:
        raise Exception("Must have balanced number of replicates per condition for now")

    n = len(df.columns) // 2
    os.system(f"mkdir -p {save_path}/tmp")

    # If exists, read results file where trial results will be concatenated to
    results, results_file, existing_trials = open_bootstrap_results(save_path, method, name)

    if results is None:
        logging.info("Initializing resultsfile")
        results_file = results_file.replace("trials*", "trials0")
        os.system(f"touch {results_file}")
        existing_trials = 0

    if existing_trials >= trials:
        logging.info(f"Already have {existing_trials} trials, returning")
        return "returned_early"

    for trial in range(existing_trials + 1, trials + 1):
        outfile_dea = f"{save_path}/tmp/tab.tmp.trial{trial}.csv"

        np.random.seed(trial)  # important for multiprocessing

        # maxiter attempts if DEA fails for small N (matrix not full rank error if too many covariates)
        for a in range(1, maxiter):
            np.random.seed(trial + (a - 1) * 1000)  # for first iteration, use trial number as seed
            # preserve matched samples
            if design == "paired":
                ind = np.array(np.random.choice(range(0, n), n))
                ind = np.concatenate([ind, ind + n])
                ind = np.sort(ind)
                df_bag = df.iloc[:, ind]

            else:
                bootstrap_samples_n = np.random.choice(df.columns[:n], n)
                bootstrap_samples_t = np.random.choice(df.columns[n:], n)
                bs = list(bootstrap_samples_n) + list(bootstrap_samples_t)
                df_bag = df[bs]

            logging.info(f"Running trial: {trial}, samples: {df_bag.columns}, path: {save_path}")

            if design == "custom" and meta is not None:
                meta_sub = meta.loc[df_bag.columns]
                meta_sub.copy()
                # add suffix to filename avoid multiprocess conflict
                design_sub = f"{save_path}/tmp/design.trial{trial}.csv"
                meta_sub.index = pd.Index([col + str(i) for i, col in enumerate(meta_sub.index)])
                meta_sub.to_csv(design_sub)
            elif design in ["paired", "unpaired"]:
                design_sub = design

            df_bag.columns = [col + str(i) for i, col in enumerate(df_bag.columns)]

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
            results = pd.concat([results, trial_results])

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

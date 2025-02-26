import sys
import os
import glob
import logging
from pathlib import Path

import pandas as pd

def main(savepath, name, trials, clean_up):
    old_merged_file = Path(f"{savepath}/{name}_trials_merged_*.csv")
    final_output = Path(f"{savepath}/{name}_trials_merged_{trials}.csv")
    matched_files = glob.glob(str(old_merged_file))
    if matched_files:
        old_merged_file = matched_files[0]
        existing_trials = int(old_merged_file.split("_trials_merged_")[1].split(".")[0])
        logger.info(f"Found: {existing_trials} existring trials, appending new ones...")
    else:
        old_merged_file = None
        logger.info("No merged file found, initializing...")
        os.system(f"touch {final_output}")

    trial_files = Path(f"{savepath}/{name}_trial_*.csv")
    matched_files = sorted(glob.glob(str(trial_files)))

    if not matched_files:
        raise Exception("No trials fond...")
    
    for tf in matched_files:
        trial = int(tf.split("_trial_")[-1].split(".csv")[0])
        logger.info(trial)
        if trial == 1:
            os.system(f"head -n 1 {tf} > {old_merged_file}") # Write header from the first trial
        os.system(f"tail -n +2 -q {tf} >> {old_merged_file}") # Append all trials, skipping headers

        if clean_up:
            os.system(f"rm {tf}")

    os.system(f"mv {old_merged_file} {final_output}")


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    savepath = sys.argv[1]
    name = sys.argv[2]
    trials = sys.argv[3]
    clean_up = sys.argv[4]

    main(savepath, name, trials, clean_up)
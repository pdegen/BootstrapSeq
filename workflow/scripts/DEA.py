import logging
import os
from pathlib import Path

import pandas as pd
import rpy2.robjects as ro
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

from R_wrappers import pd_to_R


def run_dea(df, outfile, method, overwrite=False, design="paired", lfc=0, verbose=False, **kwargs):
    """Wrapper to call appropriate R method to run differential expression analysis

    Parameters
    ----------
    df : pd.DataFrame, count data with m rows, n columns
    method: str, "edgerqlf", "edgerlrt" or "deseq2"
    overwrite: bool, overwrite existing results table if it already exists
    design: str, only use "paired" design matrix for this project
    lfc: float, formal log2 fold change threshold when testing for differential expression
    kwargs: additional keyword arguments passed to R method
    """

    script_dir = os.path.dirname(os.path.abspath(__file__))  # Get current script directory
    r_script_path = os.path.join(script_dir, "R_functions.r")  # Construct full path
    ro.r["source"](r_script_path)  # Loading the R script

    # Converting pd to R dataframe
    df_r = df if isinstance(df, ro.vectors.DataFrame) else pd_to_R(df)

    if not verbose:
        rpy2_logger.setLevel(logging.ERROR)

    if method.lower() in ["edgerqlf", "edgerlrt", "edger"]:
        logging.info(f"\nCalling edgeR in R with kwargs:\n{kwargs}\n")
        edgeR = ro.globalenv["run_edgeR"]  # Finding the R function in the script
        edgeR(df_r, str(outfile), design, overwrite, lfc=lfc, **kwargs)

    elif method.lower() == "deseq2":
        logging.info(f"\nCalling DESeq2 in R with kwargs:\n{kwargs}\n")
        if isinstance(df, ro.vectors.DataFrame):
            raise Exception("Not yet implemented for DESeq2: calling directly with df_r")
        DESeq2 = ro.globalenv["run_deseq2"]
        DESeq2(df_r, str(outfile), design, overwrite=overwrite, lfc=lfc, **kwargs)
    else:
        raise Exception(f"Method {method} not implemented")


def normalize_counts(df):
    """Use DESeq2 estimateSizeFactors to normalize a count matrix"""

    script_dir = os.path.dirname(os.path.abspath(__file__))  # Get current script directory
    r_script_path = os.path.join(script_dir, "R_functions.r")  # Construct full path
    ro.r["source"](r_script_path)  # Loading the R script

    df_r = df if isinstance(df, ro.vectors.DataFrame) else pd_to_R(df)  # Converting to R dataframe
    DESeq2 = ro.globalenv["run_deseq2"]
    sizeFactors = DESeq2(
        df_r, outfile="", design="paired", overwrite=True, print_summary=False, cols_to_keep="", size_factors_only=True
    )
    return df / list(sizeFactors)

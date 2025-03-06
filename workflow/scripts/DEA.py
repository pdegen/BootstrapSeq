import logging
import os

import pandas as pd
import rpy2.robjects as ro
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger

from R_wrappers import pd_to_r


def run_dea(
    df: pd.DataFrame,
    outfile: str,
    method: str,
    overwrite: bool = False,
    design: str = "paired",
    lfc: float = 0,
    verbose: bool = False,
    **kwargs,
) -> None:
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
    df_r = df if isinstance(df, ro.vectors.DataFrame) else pd_to_r(df)

    if not verbose:
        rpy2_logger.setLevel(logging.ERROR)

    if method.lower() in ["edgerqlf", "edgerlrt", "edger"]:
        logging.info(f"\nCalling edgeR in R with kwargs:\n{kwargs}\n")
        edger = ro.globalenv["run_edgeR"]  # Finding the R function in the script
        edger(df_r, str(outfile), design, overwrite, lfc=lfc, **kwargs)

    elif method.lower() == "deseq2":
        logging.info(f"\nCalling DESeq2 in R with kwargs:\n{kwargs}\n")
        if isinstance(df, ro.vectors.DataFrame):
            raise Exception("Not yet implemented for DESeq2: calling directly with df_r")
        deseq2 = ro.globalenv["run_deseq2"]
        deseq2(df_r, str(outfile), design, overwrite=overwrite, lfc=lfc, **kwargs)
    else:
        raise Exception(f"Method {method} not implemented")


def normalize_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Use DESeq2 estimateSizeFactors to normalize a count matrix"""

    script_dir = os.path.dirname(os.path.abspath(__file__))  # Get current script directory
    r_script_path = os.path.join(script_dir, "R_functions.r")  # Construct full path
    ro.r["source"](r_script_path)  # Loading the R script

    df_r = df if isinstance(df, ro.vectors.DataFrame) else pd_to_r(df)  # Converting to R dataframe
    deseq2 = ro.globalenv["run_deseq2"]
    size_factors = deseq2(
        df_r, outfile="", design="paired", overwrite=True, print_summary=False, cols_to_keep="", size_factors_only=True
    )
    return df / list(size_factors)

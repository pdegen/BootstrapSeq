import sys
import pandas as pd
import numpy as np
from pathlib import Path

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


def pd_to_R(df):
    """Convert pd dataframe to R dataframe"""
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_r = ro.conversion.py2rpy(df)
        return df_r


def R_to_pd(df_r):
    """Convert R dataframe to pd dataframe"""
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.conversion.rpy2py(df_r)
        return df


def filterByExpr_wrapper(inpath, outpath, design):
    """Filter low-expressed genes using edgeR's filterByExpr()
    
    Result will be saved as a csv file in outpath
    """
    ro.r['source'](str(Path('../scripts/R_functions.r')))  # Loading the R script
    R_filterByExpr = ro.globalenv['edgeR_filterByExpression']  # Finding the R function in the script
    R_filterByExpr(inpath, outpath, design)

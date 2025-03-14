from typing import Tuple

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


pretty_met = {
    "mcc": "MCC",
    "rec": "Recall",
    "rep": "Replicability",
    "prec": "Precision",
    "deg": "#DEGs",
    "terms": "#Terms",
    "std": "std",
    "kurt": "Kurtosis",
    "Spear": "Spearman correlation",
    "KL": "KL Divergence",
    "Rep": "Replicability",
    "Prec": "Precision",
    "Rec": "Recall",
}


def predict_metrics(observed_spearman: float) -> dict:
    metrics = ["Prec", "Rec", "Rep"]

    x1 = "Spear_Cohort_N5_median"
    x2 = "Spear_Cohort_N10_median"
    y1_suffix = "N5"
    y2_suffix = "N10"

    all_n = {5: (x1, y1_suffix), 10: (x2, y2_suffix)}

    dfm = pd.read_csv("../resources/degen_medo_results.csv", index_col=0)

    def linear(x, a, b):
        return a * x + b

    res_dict = {pretty_met[metric]: dict.fromkeys(all_n) for metric in metrics}
    for metric in metrics:
        for n in all_n:
            xx = all_n[n][0]
            yy = f"{metric}_{all_n[n][1]}"
            x = dfm[xx].dropna()
            y = dfm[yy].dropna()
            common = x.index.intersection(list(y.index))
            x, y = x.loc[common], y.loc[common]
            params, _ = curve_fit(linear, x, y)
            intersect = linear(observed_spearman, params[0], params[1])
            res_dict[pretty_met[metric]][n] = max(min(intersect, 1), 0)
            # print(f"{pretty_met[metric]} N{N}: {intersect:.2f}")
    return res_dict


def print_metrics(
    tab_truth, tab, fdr=0.05, return_metrics=False, return_classes=False
) -> Tuple[float, float, float] | Tuple[pd.Series, pd.Series, pd.Series, pd.Series] | None:
    common = tab_truth.index.intersection(tab.index)
    true = tab_truth.loc[common]["FDR"] < fdr
    pred = tab.loc[common]["FDR"] < fdr

    tp = true & pred
    fp = ~true & pred
    tn = ~true & ~pred
    fn = true & ~pred

    ltp, lfp, ltn, lfn = tp.sum(), fp.sum(), tn.sum(), fn.sum()

    assert ltp + lfp + ltn + lfn == len(common)

    squared = float((ltp + lfp) * (ltp + lfn) * (ltn + lfp) * (ltn + lfn))
    mcc = (ltp * ltn - lfp * lfn) / (np.sqrt(squared)) if squared else np.nan
    prec = ltp / (ltp + lfp) if ltp + lfp else np.nan
    rec = ltp / (ltp + lfn)

    print(f"MCC: {mcc:>10.2f}")
    print(f"Precision: {prec:.2f}")
    print(f"Recall: {rec:>7.2f}")
    print("===============")
    print(f"True: {true.sum():>9}")
    print(f"Pred: {pred.sum():>9}")
    print("===============")
    print(f"TP: {ltp:>11}")
    print(f"FP: {lfp:>11}")
    print(f"TN: {ltn:>11}")
    print(f"FN: {lfn:>11}")

    if return_metrics:
        return mcc, prec, rec

    if return_classes:
        return tp, fp, tn, fn

    return None

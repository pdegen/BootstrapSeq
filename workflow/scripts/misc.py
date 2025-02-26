import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


pretty_met = {"mcc": "MCC",
             "rec": "Recall",
             "rep": "Replicability",
             "prec": "Precision",
             "deg": "#DEGs",
             "terms": "#Terms",
              "std":"std",
              "kurt":"Kurtosis",
             "Spear": "Spearman correlation", "KL":"KL Divergence","Rep":"Replicability","Prec":"Precision","Rec":"Recall"
             }


def predict_metrics(observed_spearman):

    metrics = ["Prec", "Rec", "Rep"]

    x1 = "Spear_Cohort_N5_median"
    x2 = "Spear_Cohort_N10_median"
    y1_suffix = "N5"
    y2_suffix = "N10"

    all_N = {5: (x1,y1_suffix),
             10: (x2,y2_suffix)
            }
    
    dfm = pd.read_csv("../resources/degen_medo_results.csv", index_col=0)
    
    def linear(x,a,b):
        return a*x+b
    
    res_dict = {pretty_met[metric]: {N: None for N in all_N} for metric in metrics}
    for metric in metrics:
        for N in all_N:
            xx = all_N[N][0]
            yy = f"{metric}_{all_N[N][1]}"
            x = dfm[xx].dropna()
            y = dfm[yy].dropna()
            common = x.index.intersection(y.index)
            x, y = x.loc[common], y.loc[common]
            params, pcov = curve_fit(linear, x, y)
            intersect = linear(observed_spearman,params[0],params[1])
            res_dict[pretty_met[metric]][N] = intersect
            #print(f"{pretty_met[metric]} N{N}: {intersect:.2f}")
    return res_dict

def print_metrics(tab_truth, tab, FDR=0.05, return_metrics=False, return_classes=False):

    common = tab_truth.index.intersection(tab.index)
    true = tab_truth.loc[common]["FDR"]<FDR
    pred = tab.loc[common]["FDR"]<FDR
    
    TP = true & pred
    FP = ~true & pred
    TN = ~true & ~pred
    FN = true & ~pred

    lTP, lFP, lTN, lFN = TP.sum(), FP.sum(), TN.sum(), FN.sum()

    assert lTP+lFP+lTN+lFN == len(common)

    squared = float((lTP + lFP) * (lTP + lFN) * (lTN + lFP) * (lTN + lFN))
    mcc = (lTP * lTN - lFP * lFN) / (np.sqrt(squared)) if squared else np.nan
    prec = lTP / (lTP+lFP) if lTP+lFP else np.nan
    rec = lTP / (lTP+lFP)
    
    print(f"MCC: {mcc:>10.2f}")
    print(f"Precision: {prec:.2f}")
    print(f"Recall: {rec:>7.2f}")
    print("===============")
    print(f"True: {true.sum():>9}")
    print(f"Pred: {pred.sum():>9}")
    print("===============")
    print(f"TP: {lTP:>11}")
    print(f"FP: {lFP:>11}")
    print(f"TN: {lTN:>11}")
    print(f"FN: {lFN:>11}")
    
    if return_metrics:
        return mcc, prec, rec

    if return_classes:
        return TP, FP, TN, FN
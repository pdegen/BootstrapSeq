import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from scipy.optimize import curve_fit


sns.set_style("whitegrid", {"axes.linewidth": 2, "axes.edgecolor": "black"})


def make_volcano(tab, lfc=0, FDR=0.05, title="", ylim=np.inf):
    sig = tab[(tab["FDR"] < FDR) & (tab["logFC"].abs() > lfc)]
    sns.scatterplot(x=tab["logFC"], y=-np.log10(tab["FDR"]), edgecolor=None, color="grey")
    sns.scatterplot(x=sig["logFC"], y=-np.log10(sig["FDR"]), edgecolor=None)
    plt.ylabel("-log10 FDR")
    plt.axhline(-np.log10(FDR), ls="--", color="red")
    if lfc > 0:
        plt.axvline(lfc, ls="--", color="red")
        plt.axvline(-lfc, ls="--", color="red")
    if ylim < np.inf:
        plt.ylim(-0.05 * ylim, ylim)
    plt.title(f"{title} DEGs: {len(sig)}")


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
order_rep = [
    "LMAB",
    "HRLB",
    "HRLA",
    "GIPF",
    "PRAD",
    "BSHR",
    "BSLB",
    "BSLA",
    "GATB",
    "LIHC",
    "THCA",
    "BRCA",
    "KIRC",
    "LUAD",
    "COAD",
    "LUSC",
    "HSPL",
    "SNF2",
]
boxprops = dict(boxstyle="round", facecolor="#e4eaf3", alpha=1, edgecolor="#2a3b76")

palette = sns.color_palette("crest", n_colors=len(order_rep))
palette_ordered = {data: color for data, color in zip(order_rep, palette[: len(order_rep)])}

cohorts = range(50)
reference = "Cohort"
metric_suffix = "median"
metric_prefix = "Spear"  # KL, Spear

suffix = ""  # for DEGs
fit_prec = "linear"

if metric_prefix == "KL":
    fit_prec = "linear"

y_prefixes = ["Prec", "Rec", "Rep"]
x1 = f"{metric_prefix}_{reference}_N5_{metric_suffix}"
x2 = f"{metric_prefix}_{reference}_N10_{metric_suffix}"
y1_suffix = f"N5{suffix}"
y2_suffix = f"N10{suffix}"


def compare_plot(
    metric_prefix=metric_prefix,
    metric_suffix=metric_suffix,
    x1=x1,
    x2=x2,
    y1_suffix=y1_suffix,
    y2_suffix=y2_suffix,
    y_prefixes=y_prefixes,
    N1=5,
    N2=10,
    observed_spearman=None,
):
    all_N = {N1: (x1, y1_suffix), N2: (x2, y2_suffix)}

    dfm = pd.read_csv("../resources/degen_medo_results.csv", index_col=0)

    scale = 1.24
    figsize = (scale * 7.2, scale * (-1 + 4 * len(y_prefixes)))
    fig, axes = plt.subplots(len(y_prefixes), 2, figsize=figsize, sharex=False, sharey=False)

    for ax, y_prefix in zip(axes, y_prefixes):
        ax = ax.flatten()

        sns.scatterplot(
            data=dfm,
            y=f"{y_prefix}_{y1_suffix}",
            x=x1,
            hue=dfm.index,
            style=dfm.index,
            hue_order=order_rep,
            style_order=order_rep,
            s=200,
            ax=ax[0],
            palette=palette,
        )
        sns.scatterplot(
            data=dfm,
            y=f"{y_prefix}_{y2_suffix}",
            x=x2,
            hue=dfm.index,
            hue_order=order_rep,
            style_order=order_rep,
            style=dfm.index,
            s=200,
            ax=ax[1],
            palette=palette,
        )

        if fit_prec == "linear":
            sns.regplot(data=dfm, y=f"{y_prefix}_{y1_suffix}", x=x1, ax=ax[0], scatter_kws={"s": 0}, order=1)
            sns.regplot(data=dfm, y=f"{y_prefix}_{y2_suffix}", x=x2, scatter_kws={"s": 0}, ax=ax[1], order=1)

        for N, a in zip(all_N, ax):
            xx = all_N[N][0]
            yy = f"{y_prefix}_{all_N[N][1]}"
            x = dfm[xx].dropna()
            y = dfm[yy].dropna()
            common = x.index.intersection(y.index)
            x, y = x.loc[common], y.loc[common]

            if fit_prec == "binormal":

                def binormal(x, a, b):
                    return stats.norm.cdf(a * stats.norm.ppf(x) + b)

                p0 = [3, -2]  # this is an mandatory initial guess
                params, pcov = curve_fit(binormal, x, y, p0=p0, bounds=(-2.9, np.inf))
                sigma_ab = np.sqrt(np.diagonal(pcov))
                xlin = np.linspace(0.7, 1, 100)
                y_binormal = binormal(xlin, *params)
                sns.lineplot(x=xlin, y=y_binormal, color="#4d72b0", lw=2, zorder=99, ax=a)
                bound_upper = binormal(xlin, *(params + sigma_ab))
                bound_lower = binormal(xlin, *(params - sigma_ab))
                a.fill_between(xlin, bound_lower, bound_upper, color="#e4eaf3", alpha=1, zorder=0)

            else:
                r_val, p_val = stats.pearsonr(x, y)
                r2_val = r_val**2

                if y_prefix == "Prec":
                    loc = (0.95, 0.05)
                    ha = "right"
                    va = "bottom"
                else:
                    loc = (0.05, 0.95)
                    ha = "left"
                    va = "top"
                a.text(
                    loc[0],
                    loc[1],
                    f"r = {r_val:.2f}\nr² = {r2_val:.2f}\np = {p_val:.2e}",
                    transform=a.transAxes,
                    fontsize=13,
                    va=va,
                    ha=ha,
                    bbox=boxprops,
                )

            if observed_spearman and fit_prec == "linear":

                def linear(x, a, b):
                    return a * x + b

                params, pcov = curve_fit(linear, x, y)
                intersect = linear(observed_spearman, params[0], params[1])
                a.plot((observed_spearman, observed_spearman), (-1, intersect), ls="--", color="red", label="Observed")
                a.plot((-1, observed_spearman), (intersect, intersect), ls="--", color="red")

            a.set(ylabel=f"Median {pretty_met[y_prefix]}")

    ### MISC.

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles[::-1],
        labels[::-1],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        framealpha=1,
        title=None,
        ncol=1,
        markerscale=1,
    )

    for i, a in enumerate(axes.flatten()):
        a.legend().remove()
        a.set_box_aspect(1)
        a.set(xlabel=(f"{metric_suffix.split('_')[-1].capitalize()} {pretty_met[metric_prefix]}"))
        a.set(ylim=(-0.05, 1.05))
        a.set(xlim=(0.71, 0.99))
        if "Spear" in x1:
            pass
            a.xaxis.set_ticks(np.arange(0.75, 1, 0.05))
        a.annotate(
            chr(ord("A") + i),
            xy=(-0.08, 1.08),
            xycoords="axes fraction",
            weight="bold",
            va="center",
            ha="center",
            fontsize=20,
        )
        a.set_title(f"n={N1}" if i % 2 == 0 else f"n={N2}", size=16)

    return fig

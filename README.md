# BootstrapSeq

Bootstrap resample your low-powered RNA-Seq data set to estimate the expected reliability of downstream differential expression and enrichment results.

## Instructions

This pre-release version is provided as a Jupyter notebook with further instructions, see [src/bootstrapseq.ipynb](./src/bootstrapseq.ipynb). A dedicated Python package and Docker image may be released in the future (the current Dockerfile is just a placeholder).

Briefly, users provide a raw count matrix and a design matrix. Users must additionally have either edgeR or DESeq2 installed. The notebook will run bootstrapped differential expression analyses and compute the Spearman rank correlations for logFC estimates obtained from the bootstrapped and original data sets.

In Degen and Medo (2025), we show that data sets with high (>0.85) Spearman correlation have overall higher precision, recall, and replicability. The figure below shows our results for 18 different data sets.

![Fig. 5 from Degen and Medo 2025](./assets/Fig5.png)

### Number of bootstrap trials

In our original study, we limited the bootstrapping to $k=25$ trials because of the large (1'800) number of cohorts we studied. However, in real world scenarios where practitioners have a handful of data sets at best, the number of trials can be readily increased.

In a data set with $n$ biological replicates, the number of distinct bootstrap trials is ${2n-1 \choose n}$, as shown in the figure below for $n$ up to 10. For data sets with $n=5$, exhausting all possible combinations is easily achievable and should only take a few minutes when running on multiple cores.

![Combinations vs replicates](./assets/trials.png)

## To do

- Run bootstrap trials in parallel

- Let user call own R or Python script to perform DEA

- Recommend number of trials based on variability of Spearman metric

- Option to exhaust all combinations instead of random sampling

- Docker image

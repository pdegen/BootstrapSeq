# BootstrapSeq

Bootstrap resample your low-powered RNA-Seq data set to estimate the expected reliability of downstream differential expression and enrichment results. Briefly, users provide a raw count matrix and a design matrix. The provided Snakemake workflow (or Juypter notebook) will run bootstrapped differential expression analyses and compute the Spearman rank correlations for logFC estimates obtained from the bootstrapped and original data sets.

In Degen and Medo (2025), we show that data sets with a high (>0.9) Spearman correlation have overall higher precision, recall, and replicability. Conversely, data sets with a low (<0.8) correlation are prone to false positives and low replicability. The figure below shows our results for 18 different data sets.

![Fig. 5 from Degen and Medo 2025](./assets/Fig5.png)

## Instructions

In a future version, it will be possible to do the bootstrapping with a user-provided Python or R script for custom log fold change estimation. For now, edgeR is used for this. To run the workflow, we provide three options, all of which start by cloning the repository:

- `git clone https://github.com/pdegen/BootstrapSeq.git`

### Option 1: Jupyter Notebook

A Jupyter notebook with further instructions can be found in [notebooks/bootstrapseq.ipynb](notebooks/bootstrapseq.ipynb). This option does not support parallelization for now. You can handle package installation yourself or use the conda environment from Option 2.

### Option 2: Snakemake

1. Create a conda environment using [workflow/envs/environment.yaml](workflow/envs/environment.yaml)
   - `conda env create -f workflow/envs/environment.yaml`
   - `conda activate bootstrapseq`

2. Edit [config/config.yaml](config/config.yaml) as needed or create a new config and define the filepath in [workflow/Snakefile](workflow/Snakefile)

3. From the project root, run: `snakemake --cores 4` (adjust number of cores as needed)

The workflow will create a merged table with edgeR differential expression results from all trials, as well as a json file with summary statistics from calculated Spearman corelations.

### Option 3: Containerized Snakemake

For maximum reproducibility, the Snakemake workflow can also be run with [Apptainer](https://apptainer.org/docs/admin/main/installation.html) (formerly Singularity), which has to be installed separately. Then:

1. Edit [config/config.yaml](config/config.yaml) as needed or create a new config and define the filepath in [workflow/Snakefile](workflow/Snakefile)

2. From the project root, run: `snakemake --use-singularity --use-conda  --cores 4`

This command pulls the BootstrapSeq Docker image from [DockerHub](https://hub.docker.com/repository/docker/pdegen/bootstrapseq/general) with the corresponding conda environment.

### Number of bootstrap trials

In our original study, we limited the bootstrapping to $k=25$ trials because of the large (1'800) number of cohorts we studied. However, in real world scenarios where practitioners have a handful of data sets at best, the number of trials can be readily increased.

In a data set with $n$ biological replicates, the number of distinct bootstrap trials is ${2n-1 \choose n}$, as shown in the figure below for $n$ up to 10. For data sets with $n=5$, exhausting all possible combinations is easily achievable and should only take a few minutes when running on multiple cores.

![Combinations vs replicates](./assets/trials.png)

## To do

- Let user call own R or Python script to perform DEA

- Recommend number of trials based on variability of Spearman metric

- Option to exhaust all combinations instead of random sampling

- Support unbalanced designs

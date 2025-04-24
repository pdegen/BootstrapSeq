FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2962ae1422db5253d061973c1eb98f2fd4627b9f8f67c76be6883cfb5e6817bf"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/environment.yaml
#   prefix: /conda-envs/ea58bcb213908bbfc4ee610fa00ee166
#   name: bootstrapseq
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - python=3.12.5=h2ad013b_0_cpython
#     - ipykernel=6.29.5=pyh3099207_0
#     - ipython=8.27.0=pyh707e725_0
#     - bioconductor-edger=4.0.16=r43hf17093f_1
#     - matplotlib-base=3.9.2=py312hd3ec401_1
#     - seaborn=0.13.2=hd8ed1ab_2
#     - pandas=2.2.2=py312h1d6d2e6_1
#     - pyyaml=6.0.2=py312h178313f_2
#     - rpy2=3.5.11=py312r43hc7c0aa3_3
#     - r-dplyr=1.1.4=r43h0d4f4ea_1
#     - scipy=1.14.1=py312h7d485d2_0
#     - snakemake=9.3.0=hdfd78af_0
RUN mkdir -p /conda-envs/ea58bcb213908bbfc4ee610fa00ee166
COPY workflow/envs/environment.yaml /conda-envs/ea58bcb213908bbfc4ee610fa00ee166/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/ea58bcb213908bbfc4ee610fa00ee166 --file /conda-envs/ea58bcb213908bbfc4ee610fa00ee166/environment.yaml && \
    conda clean --all -y

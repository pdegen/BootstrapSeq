FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="de0db4f1667247a76cef076a5d84458c1afd024f47d8f2f9eec9d922a03295fc"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/environment.yaml
#   prefix: /conda-envs/21e68e3d865f656e3e0a7a19aa69b5a9
#   name: bootstrapseq
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - python=3.12.5=h2ad013b_0_cpython
#     - ipykernel=6.29.5=pyh3099207_0
#     - ipython=8.27.0=pyh707e725_0
#     - matplotlib-base=3.9.2=py312hd3ec401_1
#     - seaborn=0.13.2=hd8ed1ab_2
#     - pandas=2.2.2=py312h1d6d2e6_1
#     - pyyaml=6.0.2=py312h178313f_2
#     - rpy2=3.5.11=py312r43hc7c0aa3_3
#     - scipy=1.14.1=py312h7d485d2_0
#     - snakemake=8.28.0=hdfd78af_0
RUN mkdir -p /conda-envs/21e68e3d865f656e3e0a7a19aa69b5a9
COPY workflow/envs/environment.yaml /conda-envs/21e68e3d865f656e3e0a7a19aa69b5a9/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/21e68e3d865f656e3e0a7a19aa69b5a9 --file /conda-envs/21e68e3d865f656e3e0a7a19aa69b5a9/environment.yaml && \
    conda clean --all -y

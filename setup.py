from setuptools import setup, find_packages

setup(
    name="BootstrapSeq",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["jupyter", "pandas", "seaborn"],
    extras_require={"bio": ["edgeR", "DESeq2"]},
    author="Peter Degen",
    author_email="your.email@example.com",
    description="""Bootstrap resample your small-powered RNA-Seq data set to estimate the expected reliability of
    downstream differential expression and enrichment results.""",
    url="https://github.com/pdegen/BootstrapSeq",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)

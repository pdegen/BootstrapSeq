# Base image with miniconda
FROM continuumio/miniconda3

# Set working directory
WORKDIR /src

# Copy environment file and install dependencies
COPY environment.yaml .
RUN conda env create -f environment.yaml

# Activate the environment
SHELL ["conda", "run", "-n", "bootstrapseq", "/bin/bash", "-c"]

# Copy the project files
COPY src /src

# Set the default command (change main.py to your project’s main script)
CMD ["python", "main.py"]
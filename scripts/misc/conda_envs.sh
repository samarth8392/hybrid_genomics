#!/bin/bash

# Load miniconda
module load miniconda3/24.1.2-py310

# Configure conda (first time only - run these manually if needed)
# conda config --remove channels defaults
# conda config --add channels conda-forge
# conda config --set channel_priority strict

# Create environment
conda create -n bcftools_env

# Activate environment
source activate bcftools_env

# Install bcftools
conda install -c bioconda bcftools

# Verify installation
bcftools --version

# # For future use:
# module load miniconda3/24.1.2-py310
# source activate bcftools_env
# bcftools --version
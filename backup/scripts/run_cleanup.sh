#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1

# This script should be run from the results folder, using the metabar_2021 environment

Rscript "../scripts/clean_up_ASV_wise.R"
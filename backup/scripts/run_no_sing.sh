#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 1

# This script should be run from the results folder

Rscript "../scripts/no_sing_ASV_wise.R"
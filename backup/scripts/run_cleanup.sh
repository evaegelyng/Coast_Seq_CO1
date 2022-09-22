#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 1

Rscript "../scripts/clean_up_ASV_wise.R"
#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 128G
#SBATCH -c 1
#SBATCH -t 12:00:00

Rscript scripts/hmsc_COI_rich_water_plus_test.R
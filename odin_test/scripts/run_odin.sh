#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 16G
#SBATCH -c 6

# This script should be run from the results folder, using the mjolnir environment

Rscript "../scripts/odin_test.r"
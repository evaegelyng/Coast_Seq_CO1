args = commandArgs(trailingOnly=TRUE)

# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Now you can enter the total number of cores available in the system, for full computing power.
cores <- 3 # EES: Seems to be enough with a few cores

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[1])

# HELA: Hierarchical Elimination of Lurking Artifacts

# This function uses the uchime_denovo algorithm implemented in VSEARCH to remove chimaeric sequences from the dataset.
# HELA works in a sample-by-sample basis. HELA will process all individual fasta files in the current folder matching the pattern XXXX_sample_XXX.fasta.
# This allows for parallel computing, significantly decreasing calculation times.  
# The final dataset output is in VSEARCH format, so it can be directly fed into SWARM (ODIN).

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib,cores)
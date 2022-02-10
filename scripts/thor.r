args = commandArgs(trailingOnly=TRUE)

# Input path to the OTU table for the final combined library
lib_in <- as.character(args[1])

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[2])

cores <- 50
  
mjolnir5_THOR(lib,cores,tax_dir="~/taxo",ref_db="DUFA_COLR_20210723.fasta",taxo_db="taxo_NCBI_20210720")


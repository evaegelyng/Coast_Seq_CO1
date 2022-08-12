args = commandArgs(trailingOnly=TRUE)

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[1])

print(lib)

cores <- 16

# Reformat an OBITools fasta file to a vsearch file

mjolnir3_HELA_vsearch <- function(lib,cores,obipath=""){
  message("HELA will change the format to vsearch, so ODIN can use it for SWARM.")
  system(paste0("obiannotate --delete-tag=merged_sample --delete-tag=seq_rank ",lib,".new.fasta | sed 's/ count/;size/g' > ",lib,".vsearch.fasta"),intern=T,wait=T)
  message("File ",lib,".vsearch.fasta written.")
}

mjolnir3_HELA_vsearch(lib,cores)  

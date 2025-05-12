# Script to create phyloseq objects from COSQ data
# Should be run from the results folder, using the metabar_2021 environment

## Load packages
library(phyloseq)
library(dplyr)
library(stringr)

## Loading final taxonomy table 
COSQ_all <- read.table("Supplementary_table_A.txt", sep="\t", header=T, check.names=F)
nrow(COSQ_all) # 552
### Remove terrestrial taxa
COSQ_marine <- COSQ_all[which(COSQ_all$marine=="Yes" | COSQ_all$marine=="Yes (mainly Freshwater)"),]
nrow(COSQ_marine) # 484
### Remove MOTUs not identified to species level
COSQ_species <- COSQ_marine %>% filter(!str_detect(final_id_curated, 'COSQ'))
nrow(COSQ_species) # 477
### Remove MOTUs present at a single station
COSQ_prevalent <- COSQ_species[which(COSQ_species$no_stations>1),] 
nrow(COSQ_prevalent) # 368

## Write subsetted table to file
write(COSQ_prevalent$seq_id,"COSQ_pident_97_selected_250505.txt")

# Coast_Seq_CO1
Analysis of eDNA metabarcoding data from Danish coasts

After running workflow_part1.py, the following R scripts were run:

1. make_metadata.r # Produce complete metadata file
2. clean_up_OTU_wise.r # Clean the dataset based on blank controls
3. no_sing_OTU_wise.r # Remove sequences found in a single PCR replicate
4. merge_otu_table_w_classified.r # The cleaned OTU table is then merged with the taxonomy from the output file of MJOLNIR (COSQ_final_dataset)

The complete dataset was then downloaded and the remaining analysis was run locally (see "local_scripts" folder). The subset of sequences with hits of 97% similarity or more, were manually curated and normalized

To get ASV level results, OTUs that could be confidently identified to marine species, contained at least 2 ASVs and were found in at least 10 clusters were selected using select_OTUs.R. Then, continued with workflow.py, followed by:

1. asv.sh # The denoised MOTU files were combined  
2. Script2.1_DnoisE_230310.r # Further filter OTUs and remove NUMTs
3. no_sing_selected.R # Remove sequences found in a single PCR replicate
4. normalize_ASV.r # Rarefy reads to normalize sequencing depth across samples
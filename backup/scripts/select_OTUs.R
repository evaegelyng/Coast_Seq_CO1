# Script to create phyloseq objects from COSQ data
# Should be run from the results folder, using the metabar_2021 environment

## Load packages
library(phyloseq)
library(tibble)
library(dplyr)
library(ggplot2)
library(stringr)

## Loading metadata:
samples_df <- read.table("metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T, row.names=1)
## Construct phyloseq metadata table
samples=sample_data(samples_df)

## Loading final table incl. taxonomy and OTU table
COSQ_all <- read.table("NEW_pident97_data.txt", sep="\t", header=T, row.names=1,check.names=F)
### Remove sequences with no final ID, and terrestrial taxa
COSQ_marine<-COSQ_all[which(COSQ_all$final.id!="NA" & COSQ_all$Marine=="Yes"),]
### Remove sequences not identified to species level
COSQ_species<-COSQ_marine %>% filter(!str_detect(final.id, 'COSQ'))
### Count number of ASVs
sum(COSQ_species$cluster_weight)
### Count number of reads
sum(COSQ_species$total_reads)
### Subset further to MOTUS with at least 2 ASVs
COSQ<-COSQ_species[which(as.numeric(COSQ_species$cluster_weight)>=2),]
### Count number of ASVs
sum(COSQ$cluster_weight)
### Count number of reads
sum(COSQ$total_reads)

## Create OTU table
### First check where the first sample column is
COSQ[1,1:35]
### Check that the last column is a sample column
n<-ncol(COSQ)
COSQ[1,(n-1):n]
### Extract all sample columns
COSQ_otu <- COSQ[,35:n] 
### Transform to a matrix
COSQ_otu_m <- as.matrix(COSQ_otu) 
### Construct phyloseq OTU table
OTU = otu_table(COSQ_otu_m,taxa_are_rows=TRUE) 

## Create Taxonomy table
### Extract relevant columns from the COSQ table. Notes: Cols8-12 = taxonomy, Col27 = final.id
COSQ_tax <- COSQ[,c(6:12,27)] 
### Transform to a matrix
COSQ_tax_m <- as.matrix(COSQ_tax) 
### Construct phyloseq taxonomy table
TAX = tax_table(COSQ_tax_m)

## Combine metadata, OTU sample and taxonomy into one experiment-level phyloseq object
COSQ_final <- phyloseq(OTU,TAX,samples) 
COSQ_final

## Merge samples from the same sample cluster to allow estimating frequency of taxa (in terms of presence in X no. of clusters)
merged = merge_samples(COSQ_final, "cluster")
SD = merge_samples(sample_data(COSQ_final), "cluster")
print(SD)
print(merged)
sample_names(COSQ_final)
sample_names(merged)
identical(SD, sample_data(merged))
## The OTU abundances of merged samples are summed
## Check by looking at just the top10 most abundant OTUs
OTUnames10 = names(sort(taxa_sums(COSQ_final), TRUE)[1:10])
GP10 = prune_taxa(OTUnames10, COSQ_final)
mGP10 = prune_taxa(OTUnames10, merged)
C11_samples = sample_names(subset(sample_data(COSQ_final), cluster=="11"))
print(C11_samples)
otu_table(GP10)[, C11_samples]
rowSums(otu_table(GP10)[, C11_samples])
otu_table(mGP10)["11",]

## Count for each taxon the number of clusters where it is present
no_clusters<-colSums(merged@otu_table != 0)
## Count how many taxa are present in at least 5 or 10 clusters
five_c<-no_clusters[no_clusters>=5]
length(five_c)
ten_c<-no_clusters[no_clusters>=10]
length(ten_c)

## Make list of MOTUs present in at least five clusters
selection<-names(five_c)
## Extract these MOTUs from the COSQ table
COSQ_select<-COSQ[which(row.names(COSQ) %in% selection),]
## Check that the right no of MOTUs were extracted
length(selection)
nrow(COSQ_select)
### Count number of ASVs
sum(COSQ_select$cluster_weight)
### Count number of reads
sum(COSQ_select$total_reads)
## Write subsetted table to file
COSQ_select$id<-row.names(COSQ_select)
write(COSQ_select$id,"COSQ_pident_97_selected.txt")

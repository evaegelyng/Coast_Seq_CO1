# Script to normalize phyloseq objects from COSQ data
# Should be run using the metabar_2021 environment

## Load packages
library(phyloseq)
library(tibble)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
library(vegan)

## Loading metadata:
samples_df <- read.table("results/metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", header=T, row.names=1)
## Construct phyloseq metadata table
samples=sample_data(samples_df)

## Loading final ASV table
COSQ <- read.table("results/no_sing_ASV.txt", sep="\t", header=T, row.names=1,check.names=F)
### Transform to a matrix
COSQ_otu_m <- as.matrix(COSQ) 
### Construct phyloseq OTU table
OTU = otu_table(COSQ_otu_m,taxa_are_rows=TRUE) 

## Combine metadata and OTU table into one experiment-level phyloseq object
COSQ_c2 <- phyloseq(OTU,samples) 
COSQ_c2

#Remove cluster 2 (which was only sampled in 1 season)
COSQ_no_c2<-subset_samples(COSQ_c2,!cluster==2)
#Remove ASVs that are no longer represented in any samples
COSQ_no_c2 = filter_taxa(COSQ_no_c2, function(x) sum(x) > 0, TRUE)
COSQ_no_c2

# Remove OTUs with only 1 ASV
# Extract ASV table
asv_m<-data.frame(otu_table(COSQ_no_c2),check.names=F)

#Extract sample info
meta<-data.frame(sample_data(COSQ_no_c2), check.names=F)

#Reattach MOTU id for each ASV
motu<-read.table("results/COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
row.names(motu)<-motu$id
asv_motu<-merge(asv_m,motu[,c("motu","id")],by="row.names")
row.names(asv_motu)<-asv_motu$Row.names

## Remove ASV id column
asv_motu=within(asv_motu,rm("id","Row.names"))

## Count number of ASVs per MOTU
asv_count<-asv_motu%>%count(motu)
## Make a list of MOTUs with only one ASV
one_asv<-asv_count[which(asv_count$n==1),]
## Count how many MOTUs have only one ASV
nrow(one_asv)
## Remove MOTUs with only one ASV
asv_min2<-asv_motu[!asv_motu$motu %in% one_asv$motu,]
## Remove motu column
asv_min2<-within(asv_min2,rm("motu"))

## Transform to a matrix
asv_m <- as.matrix(asv_min2) 

### Construct phyloseq OTU table and metadata
ASV = otu_table(asv_m,taxa_are_rows=TRUE) 
samples = sample_data(meta)

## Combine metadata and OTU table into one experiment-level phyloseq object
COSQ_final <- phyloseq(ASV,samples) 

### Remove the PCR replicates with lowest read depth
quantile(sample_sums(COSQ_final),probs=c(0.05,.1,0.25,0.5))
COSQ_final<-prune_samples(sample_sums(COSQ_final)>2997, COSQ_final)

#Remove ASVs that are no longer represented in any samples
COSQ_final = filter_taxa(COSQ_final, function(x) sum(x) > 0, TRUE)
COSQ_final

## Rarefy PCR replicates to minimum depth

#Check read depths per PCR replicate before rarefaction
min<-min(sample_sums(COSQ_final))
min
max(sample_sums(COSQ_final))

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_merge <- rarefy_even_depth(COSQ_final, sample.size=as.numeric(min), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

### Remove ASVs from phyloseq object, that are no longer represented in any samples.
COSQ_rare = filter_taxa(COSQ_merge, function(x) sum(x) > 0, TRUE)
COSQ_rare

#Check read depths per PCR replicate after rarefaction
min(sample_sums(COSQ_rare))
max(sample_sums(COSQ_rare))

## Rarefy samples to median read depth
### First, merge PCR replicates from the same field sample
merged = merge_samples(COSQ_rare, "root")
merged

## Rebuild sample data, as the merge_samples function only handles merging of the OTU table
d<-data.frame(sample_data(merged)[,c("root","cluster","season","habitat","substrate_type","field_replicate")])

d$po<- sapply(strsplit(as.character(rownames(d)), "2C"), tail, 1)
d$pn<-gsub('\\d','', d$po)
d$pn1<-gsub(".*C(.+).*", "\\1", d$pn)
d$habitat<-ifelse(d$pn1=="EW"|d$pn1=="EB", "eelgrass", ifelse(d$pn1=="RW"|d$pn1=="RB", "rocks", "sand"))
d$substrate_type<-ifelse(grepl("B", d$pn1, fixed=T), "sediment", "water")
d$season<-ifelse(grepl("2C", as.character(rownames(d)), fixed=T), "autumn", "spring")
d$root<-rownames(d)
d$pn<-gsub('\\D','_', d$root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)
d$cluster<-as.integer(d$pn2)

sample_data(merged)<-d[,c("root","cluster","season","habitat","substrate_type","field_replicate")]

#Check read depths per sample
min(sample_sums(merged))
max(sample_sums(merged))

### Remove the samples with lowest read depth
quantile(sample_sums(merged),probs=c(0.05,0.1,0.25,0.5))
merged<-prune_samples(sample_sums(merged)>8994, merged)

#Check read depths per sample before rarefaction
min(sample_sums(merged))
max(sample_sums(merged))

### Remove ASVs from phyloseq object, that are no longer represented in any samples
COSQ_rare2 = filter_taxa(merged, function(x) sum(x) > 0, TRUE)

### Check how many ASVs are left
COSQ_rare2

## Save final file
saveRDS(COSQ_rare2,"results/COSQ_rare2_ASV_min_thres.rds")

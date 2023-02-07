# Script to normalize phyloseq objects from COSQ data
# Should be run from the results folder, using the metabar_2021 environment

## Load packages
library(phyloseq)
library(tibble)
library(dplyr)
library(ggplot2)
library(stringr)
library(vegan)

## Loading metadata:
samples_df <- read.table("metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T, row.names=1)
## Construct phyloseq metadata table
samples=sample_data(samples_df)

## Loading final table incl. taxonomy and OTU table
COSQ <- read.table("pident97_data1.txt", sep="\t", header=T, row.names=1,check.names=F)

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
### Extract relevant columns from the COSQ table. Notes: Cols8-14 = taxonomy, Col26 = score.id
COSQ_tax <- COSQ[,c(6:12,27)] 
### Transform to a matrix
COSQ_tax_m <- as.matrix(COSQ_tax) 
### Construct phyloseq taxonomy table
TAX = tax_table(COSQ_tax_m)

## Combine metadata, OTU sample and taxonomy into one experiment-level phyloseq object
COSQ_final <- phyloseq(OTU,TAX,samples) 
COSQ_final

## Rarefy PCR replicates to median depth, keeping replicates with lower depth
### Remove PCR replicates with zero reads
COSQ_final<-prune_samples(sample_sums(COSQ_final)>0, COSQ_final)

### Make a table with a column indicating which PCR replicates have a read depth above the median
readsi<-sample_sums(COSQ_final)
combinedi<-cbind(readsi, sample_data(COSQ_final))
combinedi<-data.frame(combinedi)
threshold<-round(median(combinedi$readsi))
combinedi$q<-combinedi$readsi>threshold

### Transfer the column generated above to the phyloseq object
sample_data(COSQ_final)$over_median<-combinedi$q[match(sample_data(COSQ_final)$sample_ID, combinedi$sample_ID)]

### Extract and then rarefy the PCR replicates with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(COSQ_final, over_median==TRUE), sample.size=as.numeric(threshold), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

### Extract the PCR replicates with a read depth at or below the median
below_t<-subset_samples(COSQ_final, over_median==FALSE)

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_rare<-merge_phyloseq(above_t, below_t)


## Rarefy samples to median read depth
### First, merge PCR replicates from the same field sample
merged = merge_samples(COSQ_rare, "root")

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

### Make a table with a column indicating which samples have a read depth above the median
reads<-sample_sums(merged)
combined<-cbind(reads, sample_data(merged))
combined<-data.frame(combined)
thres<-round(median(combined$reads))
combined$q<-combined$reads>thres

### Transfer the column generated above to the phyloseq object
sample_data(merged)$over_median<-combined$q[match(sample_data(merged)$root, combined$root)]

### Extract and then rarefy the samples with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(merged, over_median==TRUE), sample.size=as.numeric(thres), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

### Extract the samples with a read depth at or below the median
below_t<-subset_samples(merged, over_median==FALSE)

### Merge the rarefied samples with the low-depth samples
COSQ_rare2<-merge_phyloseq(above_t, below_t)

### Check how many OTUs are left
COSQ_rare2

## Save final files
tax_m<-data.frame(tax_table(COSQ_rare2))
otu_m<-data.frame(otu_table(COSQ_rare2),check.names=F)

write.table(data.frame(sample_data(COSQ_rare2), check.names=F), "metadata/metadata_rarefy.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(otu_m, "otu_rarefy.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "tax_rarefy.txt", sep="\t", quote=FALSE, row.names=TRUE)

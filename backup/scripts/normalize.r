# Script to normalize phyloseq objects from COSQ data
# Should be run from the results folder, using the metabar_2021 environment

## Load packages
library(phyloseq)
library(tibble)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
library(vegan)

###Load OTU table incl. taxonomy
pident70<-read.table("COSQ_final_dataset_cleaned_pident_70.tsv",sep="\t", header=T, check.names=F)

## Find first sample column
pident70[1:2,28:30] # column 29
### Check that the last column is a sample column
n<-ncol(pident70)
pident70[1,(n-1):n] 
### Extract all sample columns
otu_table <- pident70[,29:n] 

## Count no. of ASVs before removing non-marine classes
nrow(otu_table) # 7411
## Count no. of reads before removing non-marine classes
sum(colSums(otu_table)) # 102,678,511

## Load table of marine/non-marine classes
env<-read.table("18S_and_COI_classes_environment.txt",sep="\t", header=T, check.names=F)
## Check if any classes have not been assigned to marine/non-marine
pident70_na <- pident70[!pident70$class %in% env$class,]
unique(pident70_na$class) # NA

## Remove non-marine classes
pident70_mar <- pident70[pident70$class %in% env$class[env$marine=="yes"],]

### Extract all sample columns
otu_tab_mar <- within(pident70_mar,rm(marine)) 
otu_tab_mar <- otu_tab_mar[,29:n] 

## Count no. of ASVs after removing non-marine classes
nrow(otu_tab_mar) # 4078
## Count no. of reads after removing non-marine classes
sum(colSums(otu_tab_mar)) # 89,590,823

###Convert to matrix for phyloseq
otu_mat<-as.matrix(otu_tab_mar)

###Summarize no. of reads per PCR replicate
reads<-colSums(otu_mat)
mean(reads) # 24338.72
sd(reads) # 36650.16

### Extract taxonomy columns from the OTU table incl. taxonomy
tax_tab <- pident70_mar[,c("kingdom","phylum","class","order","family","genus","species","score.id")] 
# Checking that all relevant columns were included
tax_tab[1,] 
### Transform to a matrix
tax_mat <- as.matrix(tax_tab) 

#Load metadata file, containing cluster names:
metadata <- read.table("metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt")

## Combine metadata, OTU sample and taxonomy into one experiment-level phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(metadata)

COSQ_final <- phyloseq(OTU,TAX,samples) 
COSQ_final

## Rarefy PCR replicates to median depth, keeping replicates with lower depth
### Remove PCR replicates with zero reads
COSQ_final<-prune_samples(sample_sums(COSQ_final)>0, COSQ_final) # No reps removed

### Make a table with a column indicating which PCR replicates have a read depth above the median
readsi<-sample_sums(COSQ_final)
combinedi<-cbind(readsi, sample_data(COSQ_final))
combinedi<-data.frame(combinedi)
threshold<-round(median(combinedi$readsi))
combinedi$q<-combinedi$readsi>threshold

### Make histogram of raw read counts
mui <- ddply(combinedi, .(season, substrate_type), summarise, grp.mean=mean(readsi))
ggplot(combinedi, aes(x=readsi)) +
geom_histogram(aes(fill=q), position="identity", alpha=0.6, binwidth=2500) + geom_density(alpha=0.6) + geom_vline(data=mui, aes(xintercept=grp.mean), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(substrate_type ~ season, ncol=2, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Reads histogram plot", x ="Reads", y = "Count", fill = paste("Total reads > ",threshold,sep="")) + scale_x_continuous(limits=c(0,1000000)) + scale_y_continuous(limits=c(0,50))
ggsave("reads_hist_raw.pdf")

### Transfer the column generated above to the phyloseq object
sample_data(COSQ_final)$over_median<-combinedi$q[match(sample_data(COSQ_final)$sample_ID, combinedi$sample_ID)]

### Extract and then rarefy the PCR replicates with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(COSQ_final, over_median==TRUE), sample.size=as.numeric(threshold), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

# 1129OTUs were removed because they are no longer 
# present in any sample after random subsampling

### Extract the PCR replicates with a read depth at or below the median
below_t<-subset_samples(COSQ_final, over_median==FALSE)

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_merge <-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples.
COSQ_rare = filter_taxa(COSQ_merge, function(x) sum(x) > 0, TRUE)
COSQ_rare

## Rarefy samples to median read depth
### First, merge PCR replicates from the same field sample
merged = merge_samples(COSQ_rare, "sample_root")

## Rebuild sample data, as the merge_samples function only handles merging of the OTU table
d<-data.frame(sample_data(merged)[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")])

d$po<- sapply(strsplit(as.character(rownames(d)), "2C"), tail, 1)
d$pn<-gsub('\\d','', d$po)
d$pn1<-gsub(".*C(.+).*", "\\1", d$pn)
d$habitat<-ifelse(d$pn1=="EW"|d$pn1=="EB", "eelgrass", ifelse(d$pn1=="RW"|d$pn1=="RB", "rocks", "sand"))
d$substrate_type<-ifelse(grepl("B", d$pn1, fixed=T), "sediment", "water")
d$season<-ifelse(grepl("2C", as.character(rownames(d)), fixed=T), "autumn", "spring")
d$sample_root<-rownames(d)
d$pn<-gsub('\\D','_', d$sample_root)
d$pn2<-gsub(".*_(.+)__.*", "\\1", d$pn)
d$cluster<-as.integer(d$pn2)

sample_data(merged)<-d[,c("sample_root","cluster","season","habitat","substrate_type","field_replicate")]

### Make a table with a column indicating which samples have a read depth above the median
reads<-sample_sums(merged)
combined<-cbind(reads, sample_data(merged))
combined<-data.frame(combined)
thres<-round(median(combined$reads))
combined$q<-combined$reads>thres

### Transfer the column generated above to the phyloseq object
sample_data(merged)$over_median<-combined$q[match(sample_data(merged)$sample_root, combined$sample_root)]

### Extract and then rarefy the samples with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(merged, over_median==TRUE), sample.size=as.numeric(thres), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

#1396OTUs were removed

### Extract the samples with a read depth at or below the median
below_t<-subset_samples(merged, over_median==FALSE)


### Merge the rarefied samples with the low-depth samples
COSQ_merge2<-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples
COSQ_rare2 = filter_taxa(COSQ_merge2, function(x) sum(x) > 0, TRUE)

### Check how many OTUs are left
COSQ_rare2

## Save final files
tax_m<-data.frame(tax_table(COSQ_rare2))
otu_m<-data.frame(otu_table(COSQ_rare2),check.names=F)

write.table(data.frame(sample_data(COSQ_rare2), check.names=F), "metadata/metadata_rarefy_70.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(otu_m, "otu_rarefy_70.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "tax_rarefy_70.txt", sep="\t", quote=FALSE, row.names=TRUE)
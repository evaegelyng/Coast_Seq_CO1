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

## Loading metadata:
samples_df <- read.table("metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", header=T, row.names=1)
## Construct phyloseq metadata table
samples=sample_data(samples_df)

## Loading final ASV table
COSQ <- read.table("no_sing_ASV.txt", sep="\t", header=T, row.names=1,check.names=F)
### Transform to a matrix
COSQ_otu_m <- as.matrix(COSQ) 
### Construct phyloseq OTU table
OTU = otu_table(COSQ_otu_m,taxa_are_rows=TRUE) 

## Combine metadata and OTU table into one experiment-level phyloseq object
COSQ_final <- phyloseq(OTU,samples) 
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

### Make histogram of raw read counts
mui <- ddply(combinedi, .(season, substrate_type), summarise, grp.mean=mean(readsi))
ggplot(combinedi, aes(x=readsi)) +
geom_histogram(aes(fill=q), position="identity", alpha=0.6, binwidth=2500) + geom_density(alpha=0.6) + geom_vline(data=mui, aes(xintercept=grp.mean), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(substrate_type ~ season, ncol=2, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Reads histogram plot", x ="Reads", y = "Count", fill = paste("Total reads > ",threshold,sep="")) + scale_x_continuous(limits=c(0,1000000)) + scale_y_continuous(limits=c(0,50))
ggsave("reads_hist_raw_ASV.pdf")

### Transfer the column generated above to the phyloseq object
sample_data(COSQ_final)$over_median<-combinedi$q[match(sample_data(COSQ_final)$sample_ID, combinedi$sample_ID)]

### Extract and then rarefy the PCR replicates with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(COSQ_final, over_median==TRUE), sample.size=as.numeric(threshold), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

### Extract the PCR replicates with a read depth at or below the median
below_t<-subset_samples(COSQ_final, over_median==FALSE)

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_merge <-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples.
COSQ_rare = filter_taxa(COSQ_merge, function(x) sum(x) > 0, TRUE)
COSQ_rare

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
COSQ_merge2<-merge_phyloseq(above_t, below_t)

### Remove OTUs from phyloseq object, that are no longer represented in any samples
COSQ_rare2 = filter_taxa(COSQ_merge2, function(x) sum(x) > 0, TRUE)

### Check how many ASVs are left
COSQ_rare2

## Save final files
otu_m<-data.frame(otu_table(COSQ_rare2),check.names=F)

sam_dat<-data.frame(sample_data(COSQ_rare2), check.names=F)

write.table(sam_dat, "metadata/metadata_rarefy_ASV.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(otu_m, "otu_rarefy_ASV.txt", sep="\t", quote=FALSE, row.names=TRUE)

# Merge samples from the same cluster, and keep no. of positive samples instead of read counts (see Antich et al 2022)
## First add cluster information (season is included to keep dataframe format for merging)
asv_cluster<-merge(otu_m,sam_dat[,c("cluster","season")],by="row.names")
## Remove season and Row.names
asv_cluster<-within(asv_cluster,rm("season","Row.names"))
## Aggregate samples from the same cluster, counting positive samples per ASV
asv_agg<-aggregate(. ~ cluster, asv_cluster, function(x) sum(x > 0, na.rm = TRUE))

# Make NMDS ordination
## First, set row.names to cluster
row.names(asv_agg)<-asv_agg$cluster
## Remove cluster column
asv_agg=within(asv_agg,rm("cluster"))
##Remove empty columns
asv_agg<-asv_agg[, colSums(asv_agg != 0) > 0]
## Calculate Bray-Curtis distances
#dist<-vegdist(otu_agg, method="bray", binary=FALSE)
## Perform NMDS
nmds<-metaMDS(asv_agg, distance = "bray", k = 4, try = 20, trymax = 20)
## Check that stress is below 0.1 (see https://rpubs.com/CPEL/NMDS). If not, 
## try increasing k, but preferably not above 6
nmds

goodness(nmds) # Produces a results of test statistics for goodness of fit for each point

pdf("stressplot.pdf")
stressplot(nmds) # Produces a Shepards diagram
dev.off()

# Plotting points in ordination space
pdf("nmds_ASV.pdf")
plot(nmds, "sites",cex=5)   # Produces distance 
orditorp(nmds, "sites")   # Gives points labels
dev.off()

#Reattach MOTU id for each ASV
asv_agg_t<-t(asv_agg)
motu<-read.table("COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
row.names(motu)<-motu$id
asv_motu<-merge(asv_agg_t,motu[,c("motu","id")],by="row.names")
row.names(asv_motu)<-asv_motu$Row.names
## Remove cluster column
asv_motu=within(asv_motu,rm("id","Row.names"))

#Add taxonomy to ASV table
classified<-read.table("NEW_pident97_data.txt",sep="\t", header=T, check.names=F)
asv_motu$phylum<-classified$phylum[match(asv_motu$motu,classified$id)]
asv_motu$final.id<-classified$final.id[match(asv_motu$motu,classified$id)]

# Write table to a file
write.table(asv_motu,"ASV_table_clusters.tsv", sep="\t", quote=FALSE, row.names=TRUE)
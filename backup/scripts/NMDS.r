# Script to perform NMDS ordination of COSQ data
# Should be run from the results folder, using the esv environment

## Load packages
library(phyloseq)
library(tibble)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
library(vegan)

COSQ_rare2<-readRDS("COSQ_rare2_ASV.rds")

otu_m<-data.frame(otu_table(COSQ_rare2),check.names=F)

sam_dat<-data.frame(sample_data(COSQ_rare2), check.names=F)

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

#Reattach MOTU id for each ASV
asv_agg_t<-t(asv_agg)
motu<-read.table("COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
row.names(motu)<-motu$id
asv_motu<-merge(asv_agg_t,motu[,c("motu","id")],by="row.names")
row.names(asv_motu)<-asv_motu$Row.names
## Remove cluster column
asv_motu=within(asv_motu,rm("id","Row.names"))

## Count number of ASVs per MOTU
asv_count<-asv_motu%>%count(motu)
## Make a list of MOTUs with only one ASV
one_asv<-asv_count[which(asv_count$n==1),]
## Count how many MOTUs have only one ASV
nrow(one_asv)
## Remove MOTUs with only one ASV
asv_final<-asv_motu[!asv_motu$motu %in% one_asv$motu,]
## Remove MOTU id and transpose
asv=within(asv_final,rm("motu"))
asv_t<-t(asv)

## Perform NMDS
nmds<-metaMDS(asv_t, distance = "bray", k = 4, try = 20, trymax = 20)
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

## Perform NMDS at MOTU level
motu_agg<-aggregate(. ~ motu, asv_final, sum)
row.names(motu_agg)<-motu_agg$motu
motu=within(motu_agg,rm("motu"))
motu_t<-t(motu)
## Perform NMDS
nmds_motu<-metaMDS(motu_t, distance = "bray", k = 4, try = 20, trymax = 20)
## Check that stress is below 0.1 (see https://rpubs.com/CPEL/NMDS). If not, 
## try increasing k, but preferably not above 6
nmds_motu

goodness(nmds_motu) # Produces a results of test statistics for goodness of fit for each point

pdf("stressplot_motu.pdf")
stressplot(nmds_motu) # Produces a Shepards diagram
dev.off()

# Plotting points in ordination space
pdf("nmds_motu.pdf")
plot(nmds_motu, "sites",cex=5)   # Produces distance 
orditorp(nmds_motu, "sites")   # Gives points labels
dev.off()

#Add taxonomy to ASV table
classified<-read.table("NEW_pident97_data.txt",sep="\t", header=T, check.names=F)
asv_motu$phylum<-classified$phylum[match(asv_motu$motu,classified$id)]
asv_motu$final.id<-classified$final.id[match(asv_motu$motu,classified$id)]

# Write table to a file
write.table(asv_motu,"ASV_table_clusters.tsv", sep="\t", quote=FALSE, row.names=TRUE)
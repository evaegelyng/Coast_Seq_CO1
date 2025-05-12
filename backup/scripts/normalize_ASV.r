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
samples_df <- read.table("results/metadata/no_sing_ASV_metadata.txt", sep="\t", header=T, row.names=1)
## Construct phyloseq metadata table
samples <- sample_data(samples_df)

## Loading final ASV table
COSQ <- read.table("results/no_sing_ASV.txt", sep="\t", header=T, row.names=1,check.names=F)
### Transform to a matrix
COSQ_otu_m <- as.matrix(COSQ) 
### Construct phyloseq OTU table
OTU = otu_table(COSQ_otu_m,taxa_are_rows=TRUE) 

## Combine metadata and OTU table into one experiment-level phyloseq object
COSQ_c2 <- phyloseq(OTU,samples) 
COSQ_c2

## Remove eelgrass habitat which was only sampled in some clusters
COSQ_no_eel<-subset_samples(COSQ_c2,!habitat=="EB")
COSQ_no_eel<-subset_samples(COSQ_no_eel,!habitat=="EW")

# Remove ASVs that are no longer represented in any samples
COSQ_no_eel = filter_taxa(COSQ_no_eel, function(x) sum(x) > 0, TRUE)
COSQ_no_eel

#Count total remaining reads
no.eel<-data.frame(otu_table(COSQ_no_eel), check.names=F)
sum(rowSums(no.eel))

#Count no. of remaining MOTUs
#Load reference table
otu_tab<-read.table("results/COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
no.eel$motu<-otu_tab$motu[match(row.names(no.eel),otu_tab$id)]
length(unique(no.eel$motu))

#Count no. of remaining ASVs
nrow(no.eel)

## Remove cluster 2 (which was only sampled in 1 season), and
COSQ_no_c2<-subset_samples(COSQ_no_eel,!cluster==2)

# Remove ASVs that are no longer represented in any samples
COSQ_no_c2 = filter_taxa(COSQ_no_c2, function(x) sum(x) > 0, TRUE)
COSQ_no_c2

#Count total remaining reads
no.c2<-data.frame(otu_table(COSQ_no_c2), check.names=F)
sum(rowSums(no.c2))

#Count no. of remaining MOTUs
no.c2$motu<-otu_tab$motu[match(row.names(no.c2),otu_tab$id)]
length(unique(no.c2$motu))

#Count no. of remaining ASVs
nrow(no.c2)

# Remove OTUs with only 1 ASV
# Extract ASV table
#asv_m<-data.frame(otu_table(COSQ_no_c2),check.names=F)

#Extract sample info
#meta<-data.frame(sample_data(COSQ_no_c2), check.names=F)

#Reattach MOTU id for each ASV
#motu<-read.table("results/COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
#row.names(motu)<-motu$id
#asv_motu<-merge(asv_m,motu[,c("motu","id")],by="row.names")
#row.names(asv_motu)<-asv_motu$Row.names

## Remove ASV id column
#asv_motu=within(asv_motu,rm("id","Row.names"))

## Count number of ASVs per MOTU
#asv_count<-asv_motu%>%count(motu)
## Make a list of MOTUs with only one ASV
#one_asv<-asv_count[which(asv_count$n==1),]
## Count how many MOTUs have only one ASV
#nrow(one_asv)
## Remove MOTUs with only one ASV
#asv_min2<-asv_motu[!asv_motu$motu %in% one_asv$motu,]
## Remove motu column
#asv_min2<-within(asv_min2,rm("motu"))

## Transform to a matrix
#asv_m <- as.matrix(asv_min2) 

### Construct phyloseq OTU table and metadata
#ASV = otu_table(asv_m,taxa_are_rows=TRUE) 
#samples = sample_data(meta)

## Combine metadata and OTU table into one experiment-level phyloseq object
#COSQ_final <- phyloseq(ASV,samples) 
COSQ_final <- COSQ_no_c2

#Count total remaining reads
#no.sing<-data.frame(otu_table(COSQ_final), check.names=F)
#sum(rowSums(no.sing))

#Count no. of remaining MOTUs
#no.sing$motu<-otu_tab$motu[match(row.names(no.sing),otu_tab$id)]
#length(unique(no.sing$motu))

### Remove the nearly-empty samples (lowest 5%)
quantile(sample_sums(COSQ_final),probs=c(0.005,0.01,0.02,0.05))
COSQ_final<-prune_samples(sample_sums(COSQ_final)>80, COSQ_final)

#Remove ASVs that are no longer represented in any samples
COSQ_final = filter_taxa(COSQ_final, function(x) sum(x) > 0, TRUE)
COSQ_final # 3576 ASVs

#Count total remaining reads
no.low<-data.frame(otu_table(COSQ_final), check.names=F)
sum(rowSums(no.low))

#Count no. of remaining MOTUs
no.low$motu<-otu_tab$motu[match(row.names(no.low),otu_tab$id)]
length(unique(no.low$motu))

## Rarefy PCR replicates to median depth, keeping replicates with lower depth

#Check read depths per PCR replicate before rarefaction
min(sample_sums(COSQ_final))
max(sample_sums(COSQ_final))

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
ggsave("results/reads_hist_raw_ASV.pdf")

### Transfer the column generated above to the phyloseq object
sample_data(COSQ_final)$over_median<-combinedi$q[match(sample_data(COSQ_final)$sample_ID, combinedi$sample_ID)]

### Extract and then rarefy the PCR replicates with a read depth above the median
above_t<-rarefy_even_depth(subset_samples(COSQ_final, over_median==TRUE), sample.size=as.numeric(threshold), replace=FALSE, trimOTUs = TRUE, rngseed= 13072021)

### Extract the PCR replicates with a read depth at or below the median
below_t<-subset_samples(COSQ_final, over_median==FALSE)

### Merge the rarefied PCR replicates with the low-depth PCR replicates
COSQ_merge <-merge_phyloseq(above_t, below_t)

### Remove ASVs from phyloseq object, that are no longer represented in any samples.
COSQ_rare = filter_taxa(COSQ_merge, function(x) sum(x) > 0, TRUE)
COSQ_rare

#Check read depths per PCR replicate after rarefaction
min(sample_sums(COSQ_rare))
max(sample_sums(COSQ_rare))

#Count total remaining reads
rare.tab<-data.frame(otu_table(COSQ_rare), check.names=F)
sum(rowSums(rare.tab))

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

### Remove the nearly-empty samples (lowest 5%)
quantile(sample_sums(merged),probs=c(0.005,0.01,0.02,0.05))
merged<-prune_samples(sample_sums(merged)>597, merged)
merged

#Check read depths per sample before rarefaction
min(sample_sums(merged))
max(sample_sums(merged))

#Count total remaining reads
no.low2<-data.frame(otu_table(merged), check.names=F)
sum(rowSums(no.low2))

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

### Remove ASVs from phyloseq object, that are no longer represented in any samples
COSQ_rare2 = filter_taxa(COSQ_merge2, function(x) sum(x) > 0, TRUE)

### Check how many ASVs are left
COSQ_rare2

#Check read depths per sample after rarefaction
min(sample_sums(COSQ_rare2))
max(sample_sums(COSQ_rare2))

#Count total remaining reads
rare.tab2<-data.frame(otu_table(COSQ_rare2), check.names=F)
sum(rowSums(rare.tab2))

#Count no. of remaining MOTUs
rare.tab2.t<-as.data.frame(t(rare.tab2))
rare.tab2.t$motu<-otu_tab$motu[match(row.names(rare.tab2.t),otu_tab$id)]
length(unique(rare.tab2.t$motu))

## Save final file
saveRDS(COSQ_rare2,"results/COSQ_rare2_ASV_250512.rds")

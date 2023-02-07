#This script should be run from the results folder!

library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

#Load tables

###Make phyloseq object from raw data
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

## Combine OTU sample and taxonomy into one phyloseq object
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX)

#Load metadata
metadata<-read.table("metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)

##SITE INFO
c_s<-read.table("metadata/cluster_site.txt", sep="\t", header=T)
metadata$Location<-c_s$Site_name[match(metadata$cluster, c_s$cluster)]
metadata$cluster<-as.integer(metadata$cluster)
metadata$cl_se<-as.character(paste(metadata$cluster,metadata$season,sep="_"))
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)


##Now, after interpreting the figures, proceed to normalization in order to nivelate all field replicates at same depth

##Create final categorization before normalization
dir.create("depth_cutoff")

readsi<-sample_sums(DADAwang1)
combinedi<-cbind(readsi, sample_data(DADAwang1))
combinedi<-data.frame(combinedi)

quant<-c(.01, .025, .05, .1, 100)

beta_br<-expand.grid("Quantile"=quant, "Threshold"=NA, "Total_replicates_discarded"=NA, "Total_reads"=NA, "Total_ASVs"=NA, "Field_replicates"=NA, "Depth"=NA, "Mean_richness"=NA, "SD_richness"=NA, "Mean_bray_centroid_field_rep"=NA, "SD_bray_centroid_field_rep"=NA, stringsAsFactors=FALSE)

for (i in 1:length(quant)) {

threshold<-ifelse(quant[i]==100, round(median(combinedi$readsi)), round(quantile(combinedi$readsi, quant[i]), digits=-3))

beta_br[beta_br$Quantile==quant[i],"Threshold"]<-threshold

combinedi$q<-combinedi$readsi>threshold
zcom<-ddply(combinedi, .(root), transform, Total_overq=length(unique(sample_ID[q==T])))

#samples with all 4 PCR replicates above threshold
#r_et
zcom$Norm_rule<-ifelse(zcom$Total_overq==4, "r_et", NA)

#samples with 3 PCR replicates above threshold
#r_to_1
xcom<-ddply(zcom, .(root), transform, 
Norm_rule=ifelse(Total_overq==4, "r_et",
ifelse(length(unique(sample_ID[q==T]))==3&&length(unique(sample_ID[readsi>=threshold+threshold/3]))==3, "r_to_1", NA)))

#r_to_2
ycom<-ddply(xcom, .(root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&length(unique(sample_ID[readsi>=threshold+threshold/2]))==2,"r_to_2", NA)))

#r_to_3
vcom<-ddply(ycom, .(root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&max(readsi)>=2*threshold,"r_to_3", NA)))

#r_to_4
ucom<-ddply(vcom, .(root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.85*threshold&&length(unique(sample_ID[readsi>=1.05*threshold]))==3,"r_to_4", NA)))

#r_to_5
icom<-ddply(ucom, .(root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.85*threshold&&length(unique(sample_ID[readsi>=1.075*threshold]))==2,"r_to_5", NA)))

#r_to_6
ecom<-ddply(icom, .(root), transform, 
Norm_rule=ifelse(!is.na(Norm_rule), Norm_rule, ifelse(is.na(Norm_rule)&&length(unique(sample_ID[q==T]))==3&&min(readsi)>=0.75*threshold&&max(readsi)>=1.25*threshold,"r_to_6", NA)))

#Calculate max depth per field replicate
p_norm_sample_data_h<-ddply(ecom, .(root), transform, maxi=max(readsi))

#Define step-wise normalization
p_norm_sample_data<-ddply(p_norm_sample_data_h, .(sample_ID), transform, final_rule=
ifelse(is.na(Norm_rule), "trash", 
ifelse(Norm_rule=="r_et", "threshold", 
ifelse(Norm_rule=="r_to_1"&&q=="TRUE","threshold_plus_third", 
ifelse(Norm_rule=="r_to_1"&&q=="FALSE", "trash", ifelse(Norm_rule=="r_to_2"&&q=="TRUE"&&readsi>=threshold+threshold/2,"threshold_plus_half", ifelse(Norm_rule=="r_to_2"&&q=="TRUE"&&readsi<threshold+threshold/2, "threshold", ifelse(Norm_rule=="r_to_2"&&q=="FALSE", "trash", ifelse(Norm_rule=="r_to_3"&&q=="TRUE"&&readsi>=2*threshold,"threshold_twice", ifelse(Norm_rule=="r_to_3"&&q=="TRUE"&&readsi<2*threshold, "threshold", ifelse(Norm_rule=="r_to_3"&&q=="FALSE", "trash",
ifelse(Norm_rule=="r_to_4"&&q=="TRUE"&&readsi>=1.05*threshold,"threshold_1050",  ifelse(Norm_rule=="r_to_4"&&q=="FALSE", "keep_0850",
ifelse(Norm_rule=="r_to_5"&&q=="TRUE"&&readsi>=1.075*threshold,"threshold_1075", 
ifelse(Norm_rule=="r_to_5"&&q=="TRUE"&&readsi<1.075*threshold, "threshold", ifelse(Norm_rule=="r_to_5"&&q=="FALSE", "keep_0850",
ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi>=1.25*threshold&&readsi>=maxi,"threshold_1250",
ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi>=1.25*threshold&&readsi<maxi,"threshold", ifelse(Norm_rule=="r_to_6"&&q=="TRUE"&&readsi<1.25*threshold, "threshold", ifelse(Norm_rule=="r_to_6"&&q=="FALSE", "keep_0750",
"issue"))))))))))))))))))))

#update object's sample data
sample_data(DADAwang1)$final_rule<-p_norm_sample_data$final_rule[match(sample_data(DADAwang1)$sample_ID, p_norm_sample_data$sample_ID)]
steps<-levels(factor(p_norm_sample_data$final_rule))

#Check trash samples and save it
ssss<-subset(p_norm_sample_data, final_rule=="trash")

beta_br[beta_br$Quantile==quant[i],"Total_replicates_discarded"]<-length(levels(factor(ssss$sample_ID)))

#Normalize each category according to the associated sample sizes
steps=steps[steps!="trash"]

sizes<-c(threshold,
0.85*threshold,
1.075*threshold,
1.050*threshold,
threshold+threshold/3,
threshold+threshold/2,
threshold*2,
0.75*threshold,
1.250*threshold
)

sizes2<-round(as.numeric(sizes),digits=0)

prtc<-c("threshold",
"keep_0850",
"threshold_1075",
"threshold_1050",
"threshold_plus_third",
"threshold_plus_half",
"threshold_twice",
"keep_0750",
"threshold_1250")

prtc2<-as.character(prtc)
gaba<-data.frame(cbind(sizes2, prtc2), stringsAsFactors = FALSE)
pepe<-gaba[match(steps, gaba$prtc2),]
stnd<-subset(pepe, !is.na(pepe$prtc2))

for(e in 1:nrow(stnd)) {
  assign(paste0("samples", stnd[e,2]), rarefy_even_depth(subset_samples(DADAwang1, final_rule==stnd[e,2]), sample.size=as.numeric(stnd[e,1]), replace=FALSE, rngseed= 13072021))
}

#merge objects
myobs<-paste0("samples",stnd$prtc2)
tudao_r_et<-do.call(merge_phyloseq, mget(myobs), quote = FALSE)

tudao_pos_merge<-merge_samples(tudao_r_et, "root")

#Re-build sample_data
tudao0 = filter_taxa(tudao_pos_merge, function(x) sum(x) > 0, TRUE)
p_silva = prune_samples(sample_sums(tudao0)>0,tudao0)

otug<-otu_table(p_silva)
otugi<-t(data.frame(otug, check.names=F))
rich<-colSums(otugi != 0)

#Write info
beta_br[beta_br$Quantile==quant[i],"Total_reads"]<-sum(sample_sums(p_silva))
beta_br[beta_br$Quantile==quant[i],"Total_ASVs"]<-ntaxa(p_silva)
beta_br[beta_br$Quantile==quant[i],"Field_replicates"]<-nsamples(p_silva)
beta_br[beta_br$Quantile==quant[i],"Depth"]<-ifelse(length(levels(as.factor(sample_sums(p_silva))))==1, levels(as.factor(sample_sums(p_silva))), mean(sample_sums(p_silva)))
beta_br[beta_br$Quantile==quant[i],"Mean_richness"]<-round(mean(rich), digits=0)
beta_br[beta_br$Quantile==quant[i],"SD_richness"]<-round(sd(rich), digits=0)

d<-data.frame(sample_data(p_silva)[,c("root","cluster","season","habitat","substrate_type","field_replicate")])

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

#Calculate beta_diversity (as PERMANOVA results)
dist_n<-vegdist(data.frame(otug), method="bray")
d$s_st_h_c<-paste(d$cluster,d$season,d$substrate_type,d$habitat, sep="_")
zin<-betadisper(dist_n, as.factor(d$s_st_h_c))
tre<-with(zin, tapply(distances, as.factor(d$s_st_h_c), "mean"))

#Write beta
beta_br[beta_br$Quantile==quant[i],"Mean_bray_centroid_field_rep"]<-mean(tre)
beta_br[beta_br$Quantile==quant[i],"SD_bray_centroid_field_rep"]<-sd(tre)

sample_data(p_silva)<-d[,c("root","cluster","season","habitat","substrate_type","field_replicate")]

p_silva.2 <- transform_sample_counts(p_silva, function(x) sqrt(x))
ordb<-ordinate(p_silva.2, "PCoA", "bray")

plot<-plot_ordination(p_silva.2, ordb, type = "sample", color="habitat", shape="season") + geom_text(aes(label=cluster), size=1.5, color="black", alpha=0.9) + geom_point(size=1.6, alpha=0.2) + stat_ellipse(aes(lty=season), type = "norm") + facet_wrap(~substrate_type, ncol=1) + theme_bw()
ggsave(filename = paste("quantile",quant[i],"ord_bc_root_silva.pdf", sep="_"), path = "depth_cutoff", plot = plot)

}


#####Last round - using media and keeping size of samples with depth below threshold
yul<-DADAwang1
sample_data(yul)$over_median<-p_norm_sample_data$q[match(sample_data(yul)$sample_ID, p_norm_sample_data$sample_ID)]

above_t<-rarefy_even_depth(subset_samples(yul, over_median==TRUE), sample.size=as.numeric(round(median(combinedi$readsi))), replace=FALSE, rngseed= 13072021)
below_t<-subset_samples(yul, over_median==FALSE)

#merge objects
tudao_r_new<-merge_phyloseq(above_t, below_t)
tudao_pos_merge_new<-merge_samples(tudao_r_new, "root")

#Re-build sample_data
tudaok = filter_taxa(tudao_pos_merge_new, function(x) sum(x) > 0, TRUE)
p_silva = prune_samples(sample_sums(tudaok)>0,tudaok)

otug<-otu_table(p_silva)
otugi<-t(data.frame(otug, check.names=F))
rich<-colSums(otugi != 0)

#Write info
beta_br[6,"Quantile"]<-"keep_all"
beta_br[6,"Total_reads"]<-sum(sample_sums(p_silva))
beta_br[6,"Total_ASVs"]<-ntaxa(p_silva)
beta_br[6,"Field_replicates"]<-nsamples(p_silva)
beta_br[6,"Depth"]<-ifelse(length(levels(as.factor(sample_sums(p_silva))))==1, levels(as.factor(sample_sums(p_silva))), mean(sample_sums(p_silva)))
beta_br[6,"Mean_richness"]<-round(mean(rich), digits=0)
beta_br[6,"SD_richness"]<-round(sd(rich), digits=0)
beta_br[6,"Total_replicates_discarded"]<-length(levels(factor(combinedi$sample_ID)))-nsamples(tudao_r_new)

d<-data.frame(sample_data(p_silva)[,c("root","cluster","season","habitat","substrate_type","field_replicate")])

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


#Calculate beta_diversity (as PERMANOVA results)
dist_n<-vegdist(data.frame(otug), method="bray")
d$s_st_h_c<-paste(d$cluster,d$season,d$substrate_type,d$habitat, sep="_")
zin<-betadisper(dist_n, as.factor(d$s_st_h_c))
tre<-with(zin, tapply(distances, as.factor(d$s_st_h_c), "mean"))

beta_br[6,"Mean_bray_centroid_field_rep"]<-mean(tre)
beta_br[6,"SD_bray_centroid_field_rep"]<-sd(tre)

sample_data(p_silva)<-d[,c("root","cluster","season","habitat","substrate_type","field_replicate")]

p_ncbi.2 <- transform_sample_counts(p_silva, function(x) sqrt(x))
ordb<-ordinate(p_ncbi.2, "PCoA", "bray")

plot<-plot_ordination(p_ncbi.2, ordb, type = "sample", color="habitat", shape="season") + geom_text(aes(label=cluster), size=1.5, color="black", alpha=0.9) + geom_point(size=1.6, alpha=0.2) + stat_ellipse(aes(lty=season), type = "norm") + facet_wrap(~substrate_type, ncol=1) + theme_bw()
ggsave("keep_below_ord_bc_root_silva.pdf",path="depth_cutoff")

write.table(beta_br, "depth_cutoff/summary_depth_cutoff_silva.txt", sep="\t", quote=FALSE)

#Save final files

tax_m<-data.frame(tax_table(p_silva))
otu_m<-data.frame(otu_table(p_silva),check.names=F)

write.table(data.frame(sample_data(p_silva), check.names=F), "metadata/f_silva_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(otu_m, "f_otu_silva.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_m, "f_tax_silva.txt", sep="\t", quote=FALSE, row.names=TRUE)



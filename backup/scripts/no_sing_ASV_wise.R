library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

#This script should be run from the "results" folder. 

#Load tables
otu_mat<-as.matrix(read.table("cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

OTUt = otu_table(otu_mat, taxa_are_rows = TRUE)
p_DADAwangt = phyloseq(OTUt)

#Load metadata
metadata<-read.table("metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))
DADAwang1 = merge_phyloseq(p_DADAwangt, sampledata)

#Filter samples and taxa with 0 reads
p_dsa = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
dsa = prune_samples(sample_sums(p_dsa)>0,p_dsa)
sum(sample_sums(DADAwang1))==sum(sample_sums(dsa))
cat("\n")
cat("total_reads_before_control_reads_removal")
sum(sample_sums(dsa))
cat("\n")

#Exclude all control samples (NTCs, CNEs and field controls)
p_velux<-subset_samples(dsa,source=="Field_sample")
p_tre = filter_taxa(p_velux, function(x) sum(x) > 0, TRUE)
pert<-p_tre
new_otu_mat<-as.matrix(data.frame(otu_table(pert), check.names=F))

pert
cat("\n")
cat("total_reads_before_singleton_removal")
sum(sample_sums(pert))
cat("\n")
cat("example_sample_with_singleton_COSQ_000000089_2C10EB3_3_before_removal")
new_otu_mat["COSQ_000000089","2C10EB3_3"]
sdptre<-data.frame(sample_data(pert))

for (i in unique(sdptre$root))
{
     idxName = grepl(unlist(i), row.names(sdptre))
     if(sum(idxName)>1)
       checkValues = rowSums(1*(new_otu_mat[,idxName]>0))
     if(sum(idxName)==1)
       checkValues = 1*(new_otu_mat[,idxName]>0)
     idxUnderThres = which(checkValues < 2)
     new_otu_mat[idxUnderThres,idxName] = 0
}

cat("\n")
cat("example_sample_with_singleton_COSQ_000000089_2C10EB3_3_after_removal")
new_otu_mat["COSQ_000000089","2C10EB3_3"]

otu_table(p_tre)<-otu_table(new_otu_mat, taxa_are_rows = TRUE)
tre = filter_taxa(p_tre, function(x) sum(x) > 0, TRUE)
tsa = prune_samples(sample_sums(tre)>0,tre)

cat("\n")
cat("total_reads_after_singleton_removal")
sum(sample_sums(tsa))
cat("\n")
tsa
shablaw<-as.matrix(data.frame(otu_table(tsa), check.names=F))

#Write files

write.table(shablaw, "no_sing_cleaned_otu_table_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(data.frame(sample_data(tsa)), "metadata/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)


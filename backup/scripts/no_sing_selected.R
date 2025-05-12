# Remember to activate conda environment "esv"

library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

#Load table
otu_tab<-read.table("results/COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
#Count no. of MOTUs before filtering
length(unique(otu_tab$motu))

#Set row names to ASV ID
row.names(otu_tab)<-otu_tab$id

#Select sample columns only
##Find first sample column
otu_tab[1:2,1:5]
##Find last sample column
n<-ncol(otu_tab)
otu_tab[1:2,(n-2):n]
##Subset to sample columns
otu_table<-otu_tab[,4:(n-2)]

#Change to phyloseq object
otu_mat<-as.matrix(otu_table)

OTUt = otu_table(otu_mat, taxa_are_rows = TRUE)
p_DADAwangt = phyloseq(OTUt)

#Load metadata
metadata<-read.table("results/metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))
DADAwang1 = merge_phyloseq(p_DADAwangt, sampledata)

#Filter samples and taxa with 0 reads
p_dsa = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
dsa = prune_samples(sample_sums(p_dsa)>0,p_dsa)
sum(sample_sums(DADAwang1))==sum(sample_sums(dsa))
cat("\n")
cat("total_reads_before_singleton_removal")
sum(sample_sums(dsa))
cat("\n")

cat("example_sample_with_singleton_COSQ_000030162_2C10EB1_4_before_removal")
new_otu_mat<-as.matrix(data.frame(otu_table(dsa), check.names=F))
new_otu_mat["COSQ_000030162","2C10EB1_4"]
sdptre<-data.frame(sample_data(dsa))

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
cat("example_sample_with_singleton__COSQ_000030162_2C10EB1_4_before_removal")
new_otu_mat["COSQ_000030162","2C10EB1_4"]

otu_table(dsa)<-otu_table(new_otu_mat, taxa_are_rows = TRUE)
tre = filter_taxa(dsa, function(x) sum(x) > 0, TRUE)
tsa = prune_samples(sample_sums(tre)>0,tre)

cat("\n")
cat("total_reads_after_singleton_removal")
sum(sample_sums(tsa))
cat("\n")
tsa
shablaw<-as.matrix(data.frame(otu_table(tsa), check.names=F))

#Write files

write.table(shablaw, "results/no_sing_ASV.txt", sep="\t", quote=FALSE, row.names=TRUE)

write.table(data.frame(sample_data(tsa)), "results/metadata/no_sing_ASV_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

# Count no. of remaining MOTUs
no.sing<-data.frame(otu_table(tsa), check.names=F)
no.sing$motu<-otu_tab$motu[match(row.names(no.sing),otu_tab$id)]
length(unique(no.sing$motu))

# Count no. of remaining ASVs
nrow(no.sing)
library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")

#This script should be run from the "results" folder. 

#import data with check.names=F
mjolnir_output<-read.table("COSQ_final_dataset.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Check where sample columns begin and end, and extract these columns. 
n<-ncol(mjolnir_output)
mjolnir_output[1,1:30] 
mjolnir_output[1,(n-1):n]
otu_mat <- mjolnir_output[ -c(n,1:26) ]

###Make phyloseq object from raw data
f_otu_mat<-as.matrix(sapply(otu_mat, as.numeric))
rownames(f_otu_mat)<-rownames(otu_mat)
OTU = otu_table(f_otu_mat, taxa_are_rows = TRUE)
p_DADAwang = phyloseq(OTU)

#From sample_names generate metadata
sample_ID<-sample_names(p_DADAwang)

###Get sample source
source<-NA
metadata<-data.frame(cbind(sample_ID, source))
rownames(metadata)<-metadata$sample_ID
metadata$p_source<- sapply(strsplit(as.character(metadata$sample_ID), "_"), head, 1)
levels(factor(metadata$p_source))

#In the below, have removed "2SN" and "2WN", as these are included by grepping "SN" and "WN".
metadata$source<-ifelse(metadata$p_source=="CNE"|metadata$p_source=="2CCNE", 
as.character("CNE"), ifelse(metadata$p_source=="control", as.character("control"), 
ifelse(grepl("SN", metadata$p_source, fixed=T)|grepl("WN", metadata$p_source, fixed=T),as.character("NTC"),
as.character("Field_sample"))))

levels(factor(metadata$source))

#Get root
head(metadata)

metadata$root<-gsub("_+[^_]+$", "",metadata$sample_ID) 

levels(factor(metadata$root))

###Get sample PCR replicate and seq run
metadata$po<- sapply(strsplit(as.character(metadata$sample_ID), "_"), tail, 1)
levels(factor(metadata$po))

metadata$PCR_replicate<-substr(metadata$po,1,1)

metadata$seq_run<-substr(metadata$po,3,3)

head(metadata)

###Get sample cluster
metadata$pn<-gsub('\\D','_', metadata$sample_ID)
metadata$pn2<-gsub(".*_(.+)__.*", "\\1", metadata$pn)
metadata[,c("root","pn2")]
levels(factor(metadata$pn2))

metadata$pn3<-gsub(".*[C]([^.]+)[_].*", "\\1", metadata$sample_ID)
metadata[,c("root","pn3")]

metadata$cluster<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn2), ifelse(metadata$source=="control", as.integer(gsub('\\D','', metadata$pn3)), NA))

metadata[,c("root","cluster")]
levels(factor(metadata$cluster))
str(metadata)

###Get field replicate
metadata$pn6<- sapply(strsplit(as.character(metadata$pn), "__"), tail, 1)
metadata$pn7<-sapply(strsplit(as.character(metadata$pn6), "_"), head, 1)

metadata$field_replicate<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn7), NA)

metadata[,c("root","field_replicate")]
levels(factor(metadata$field_replicate))
str(metadata)

###Get habitat
metadata$habitat<-ifelse(metadata$source=="Field_sample", as.character(gsub('\\d','', metadata$pn3)), NA)

metadata[,c("root","habitat")]
levels(factor(metadata$habitat))
str(metadata)

###Get extraction refs
extraction_refs<-read.table("~/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/backup/data/raw_data/extraction_refs_both_seasons.txt", sep="\t", header=T, row.names=1)
extraction_refs$sample_root<-row.names(extraction_refs)
extraction_refs$extraction_refs<-sapply(strsplit(as.character(extraction_refs$extraction_number), "_"), head, 1)

metadata$extraction_refs<-extraction_refs$extraction_refs[match(metadata$root, extraction_refs$sample_root)]

metadata[,c("root","extraction_refs")]
levels(factor(metadata$extraction_refs))
head(metadata)

###Get PSU refs
PSU_refs<-read.table("metadata/COSQ_metadata_reps.tsv", sep="\t", header=T)
head(PSU_refs)

metadata$PSU_refs<-PSU_refs$PSU[match(metadata$root, PSU_refs$root)]

metadata[,c("root","PSU_refs")]
levels(factor(metadata$PSU_refs))
head(metadata)

###Get substrate_type
metadata$substrate_type<-PSU_refs$substrate_type[match(metadata$root, PSU_refs$root)]

metadata[,c("root","substrate_type")]
levels(factor(metadata$substrate_type))
str(metadata)

###Get season
metadata$season<-PSU_refs$season[match(metadata$root, PSU_refs$root)]

metadata[,c("root","season")]
levels(factor(metadata$season))
head(metadata)

###Discard tmp variables and inspect metadata sheet
colnames(metadata)
test_md<-metadata[,c("sample_ID", "root", "source", "season", "seq_run", "substrate_type", "cluster", "habitat", "field_replicate", "extraction_refs", "PSU_refs", "PCR_replicate")]
head(test_md)
colnames(test_md)[2]<-"sample_root"

####Check NAs. Only NTCs should have no extraction refs, and NTCs should have "NA" for 
####the variables habitat, cluster, season and field replicate.
na<-test_md[is.na(test_md$extraction_refs),]
levels(factor(na$source))
levels(factor(na$habitat))
levels(factor(na$cluster))
levels(factor(na$season))
levels(factor(na$field_replicate))

write.table(test_md, "metadata/COSQ_metadata_complete.txt", sep="\t", row.names = T, quote=FALSE)
write.table(f_otu_mat, "COSQ_otu_phyloseq.txt", sep="\t", quote=FALSE, row.names=TRUE)
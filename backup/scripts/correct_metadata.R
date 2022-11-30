setwd("/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/results")

#Read metadata table made from ngsfilter files
metadata<-read.table("metadata/COSQ_metadata.tsv", sep="\t", header=T)

#Read table with corrected replicate numbers and resequencing indicated
reps<-read.table("/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/backup/data/raw_data/library_names_reps.csv",sep="\t",header=T)

#Make a column of library names in metadata to allow merging with the reps dataframe
metadata$library<-sapply(strsplit(as.character(metadata$mjolnir_agnomens), "_"), head, 1)

#Merge metadata with reps by library column
meta_reps<-merge(metadata,reps,by="library")

#Saving the sample names from the original ngsfilter files (made from tags files)
meta_reps$ngsfilter_names<-meta_reps$original_samples

#Removing replicate number from old sample names, and saving root name to a new variable
meta_reps$original_samples<-gsub("_+[^_]+$", "",meta_reps$original_samples)
meta_reps$root<-gsub("_+[^_]+$", "",meta_reps$original_samples)

#Pasting corrected replicate numbers onto sample names
meta_reps$original_samples <- apply( meta_reps[ , c(3,4) ] , 1 , paste , collapse = "_" )

#Paste sequencing replicate (original sequencing = 1s, resequencing = 2s and 3s) to new sample names
meta_reps$original_samples <- apply( meta_reps[ , c(3,6) ] , 1 , paste , collapse = "" )

#Make and export new metadata table from meta_reps
metadata_new<-meta_reps[c(2,3)]
write.table(metadata_new, "COSQ_metadata_new.tsv", sep="\t", quote=FALSE, row.names=F)

#Make a column of original replicate number
meta_reps$ngsfilter_rep<- sapply(strsplit(as.character(meta_reps$ngsfilter_names), "_"), tail, 1)

#Make a column showing whether the corrected replicate number is different from the original
meta_reps$old_new_rep_diff<-ifelse(meta_reps$replicate-as.numeric(meta_reps$ngsfilter_rep)>0,1,0)

#Make a column of substrate type from library number
meta_reps$substrate_type<-ifelse(substr(meta_reps$library,1,1)=="S","sediment","water")

#Make a column of season from library number
meta_reps$season<-ifelse(substr(meta_reps$library,2,2)==1,"spring","autumn")

#Make a column of PSU number from library number
meta_reps$PSU<-paste("PSU",substr(meta_reps$library,3,3),sep="")

#Write results to a table 
write.table(meta_reps, "COSQ_metadata_reps.tsv", sep="\t", quote=FALSE, row.names=F)

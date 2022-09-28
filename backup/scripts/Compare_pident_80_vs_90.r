#Script to compare taxonomic identifications between datasets w. similarity threshold of 80% and 90%

#Import dataset w. 90% threshold
p90<-read.table("COSQ_final_dataset_cleaned_pident_90.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Subset to the columns "alternatives" and "score.id"
sub90<-p90[,c(22,26)]
head(sub90,n=2)

#Import dataset w. 80% threshold
p80<-read.table("COSQ_final_dataset_cleaned.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Subset to the columns "alternatives" and "score.id"
sub80<-p80[,c(22,26)]
head(sub80,n=2)
#Rename columns to distinguish from 90% file, before merging the files
names(sub80)[1]<-"alternatives.pident.80"
names(sub80)[2]<-"score.id.pident.80"

#Merge columns by sequence id
p80_p90<-merge(sub80,sub90,by="row.names")
head(p80_p90,n=2)
names(p80_p90)[1]<-"id"

#Make a column indicating whether there is a difference between alternatives in the 80% vs 90% file
p80_p90$diff.altern<-ifelse(p80_p90$alternatives.pident.80==p80_p90$alternatives,"no","yes")
#Make a column indicating whether there is a difference between score.id in the 80% vs 90% file
p80_p90$diff.score.id<-ifelse(p80_p90$score.id.pident.80==p80_p90$score.id,"no","yes")

#Subset to sequences (rows) that differ between the 80% and 90% files in alternatives or score.id
diffs<-p80_p90[which(p80_p90$diff.altern=="yes"|p80_p90$diff.score.id=="yes"),]
head(diffs,n=5)
#Reorder columns to make the file easier to read
diffs_reorder<-diffs[,c(1,4,2,5,3,6,7)]
head(diffs_reorder,n=2)
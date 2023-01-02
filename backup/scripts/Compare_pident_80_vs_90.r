#Script to compare taxonomic identifications between datasets w. similarity threshold of 70% and 80%

#Import dataset w. 80% threshold
p80<-read.table("COSQ_final_dataset_cleaned_pident_80.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Subset to the columns "alternatives" and "score.id"
sub80<-p80[,c(22,26)]
head(sub80,n=2)

#Import dataset w. 70% threshold
p70<-read.table("COSQ_final_dataset_cleaned.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Subset to the columns "alternatives" and "score.id"
sub70<-p70[,c(22,26)]
head(sub70,n=2)
#Rename columns to distinguish from 80% file, before merging the files
names(sub70)[1]<-"alternatives.pident.70"
names(sub70)[2]<-"score.id.pident.70"

#Merge columns by sequence id
p70_p80<-merge(sub70,sub80,by="row.names")
head(p70_p80,n=2)
names(p70_p80)[1]<-"id"

#Make a column indicating whether there is a difference between alternatives in the 70% vs 80% file
p70_p80$diff.altern<-ifelse(p70_p80$alternatives.pident.70==p70_p80$alternatives,"no","yes")
#Make a column indicating whether there is a difference between score.id in the 70% vs 80% file
p70_p80$diff.score.id<-ifelse(p70_p80$score.id.pident.70==p70_p80$score.id,"no","yes")

#Subset to sequences (rows) that differ between the 70% and 80% files in alternatives or score.id
diffs<-p70_p80[which(p70_p80$diff.altern=="yes"|p70_p80$diff.score.id=="yes"),]
head(diffs,n=5)
#Reorder columns to make the file easier to read
diffs_reorder<-diffs[,c(1,4,2,5,3,6,7)]
head(diffs_reorder,n=2)

write.table(diffs_reorder,"diffs_p70_p80.tsv",sep="\t",row.names=F)
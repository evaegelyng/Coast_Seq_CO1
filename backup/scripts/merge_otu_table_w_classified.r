#This script should be run from the "results" folder. 

#import classification data with check.names=F
mjolnir_output<-read.table("COSQ_final_dataset.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Extract the classification columns and other non-sample columns. 
n<-ncol(mjolnir_output)
tax_mat <- mjolnir_output[ c(1:26,n) ]

#Import cleaned OTU table
otu_mat<-read.table("no_sing_cleaned_otu_table_ASV_wise.txt",sep="\t",check.names=F)

#Merge the two tables
tax_otu<-merge(tax_mat,otu_mat,by="row.names")

#Rename sequence ID column
names(tax_otu)[1]<-"id"
write.table(tax_otu,"COSQ_final_dataset_cleaned_pident70.tsv",sep="\t",row.names=F)

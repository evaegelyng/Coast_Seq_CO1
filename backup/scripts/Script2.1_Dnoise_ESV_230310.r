library("optparse")
library("stringr")
library("dplyr")
library("phyloseq")
source('scripts/Script2.2_numts.R')

output_dir = "tmp/finalfiles"

print('loading data')

ESV_data_initial <- read.csv("tmp/ESV_Adcorr.csv")

# Load file with MJOLNIR agnomens and corresponding original names
metadata <- read.table("results/metadata/COSQ_metadata_new.tsv",sep="\t",head=T,stringsAsFactors = F)
#sample_metadata_sorted <- read.table(metadata2) # this would be for the metadata in an sorted way
original_data <- read.csv("results/COSQ_Curated_LULU.tsv", sep='\t')
# Not sure if we should include only metazoans as below
#motu_taxa <- data.frame('id' = original_data$id, 'Metazoa' = c(original_data$kingdom == 'Metazoa'))

print('data loaded')

#Select sample columns only
sample_cols <- grep("sample",colnames(ESV_data_initial))
#Remove the prefix "sample." from sample names
sample_names <- gsub("sample\\.","",names(ESV_data_initial[sample_cols]))

# Change agnomens to original names
new_sample_names <- metadata$original_samples[match(sample_names,metadata$mjolnir_agnomens)]
# The below did not reveal any NAs, so did not run
#new_sample_names[is.na(new_sample_names)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names_ESV)))))
names(ESV_data_initial)[sample_cols] <- new_sample_names

original_data <- original_data[,-grep('sample',colnames(original_data))]

neg_samples <- grep("CNE|SN|WN|control", colnames(ESV_data_initial))
neg_data <- ESV_data_initial[,neg_samples]

# Remove negative controls from data
ESV_data_initial <- ESV_data_initial[,-neg_samples]

# Merge data from the same sample but different sequencing runs
# First rename row names to ESV id
row.names(ESV_data_initial)<-ESV_data_initial$id

## Create OTU table
### First check where the first sample column is
ESV_data_initial[1,1:5]
### Check where the last sample column is
n<-ncol(ESV_data_initial)
ESV_data_initial[1,(n-1):n]
### Extract all sample columns
otu <- ESV_data_initial[,5:(n-1)] 
### Transpose
otu_t <- t(otu) 
otu_t <- data.frame(otu_t)
otu_t$root=factor(substr(row.names(otu_t),1,nchar(row.names(otu_t))-2))
otu_agg <- aggregate(.~root,data=otu_t,sum)
row.names(otu_agg)<-otu_agg$root
otu_agg<-within(otu_agg,rm("root"))
otu_agg_t <- t(otu_agg)
read_data <- data.frame(otu_agg_t,check.names=F)

# Merge new OTU table with taxonomy and sequences
initial_data <- ESV_data_initial[,c(1,2,4)]

seq_data <- ESV_data_initial$sequence

filtered_data <- cbind(initial_data, read_data, seq_data)

rm(ESV_data_initial, initial_data, read_data, seq_data)

#####
# print('starting Filter 1: remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance remove blank and NEG samples')
# # Filter 1. remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance
# # remove blank and NEG samples
# data_neg_filt_deleted <- filtered_data[rowSums(neg_data)/rowSums(cbind(filtered_data[,grep("2018",colnames(filtered_data))],neg_data)) > 0.1,]
# filtered_data <- filtered_data[rowSums(neg_data)/rowSums(cbind(filtered_data[,grep("2018",colnames(filtered_data))],neg_data)) <= 0.1,]
# 
# # Add samples names and write data into neg_filtrates directory
# write.csv(filtered_data, file = paste0(output_dir,"ESV_negfilt.csv"),row.names = F)
# 
# write.csv(data_neg_filt_deleted, file = paste0(output_dir,"ESV_negfilt_deleted.csv"),row.names = F)
# rm(data_neg_filt_deleted)
# 
# print('Filter 1 finished')


#####
print('starting Filter 2: Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below 0.005% of the total reads of this sample')
# Filter 2. Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below 0.005% of the total reads of this sample

modified_ESV <- list()

relabund <- function(x){
  x/sum(x)
}



for (i in sample_metadata_sorted$sample_ID) {
  a <- relabund(filtered_data[,grep(i,colnames(filtered_data))])<0.00005 
  zeros <- filtered_data[,grep(i,colnames(filtered_data))] == 0 
  
  changes <- (a+zeros == 1)
  
  filtered_data[changes,grep(i,colnames(filtered_data))] <- 0
  
  # modified motus are stored
  modified_ESV[i] <- list(as.character(filtered_data$id[changes]))
}


# write modified MOTUs info

write.csv(filtered_data, file = paste0(output_dir,"ESV_relabundfilt.csv"),row.names = F)

sink(paste0(output_dir,"ESV_relabundfilt_modified.txt"))
print(modified_ESV)
sink()

rm(modified_ESV)

print('Filter 2 finished')


#####
print('starting Filter 3: apply minimum abundance of 5 reads')
# Filter 3. apply minimum abundance of 5 reads
filtered_data$count <- rowSums(filtered_data[,grep('2018',colnames(filtered_data))])
filtered_data <- filtered_data[filtered_data$count>4,]
colnames(filtered_data)[grep('seq_data',colnames(filtered_data))] <- 'seq'

print('Filter 3 finished')

#####
print('starting Filter 4: remove sequences of bad length')
# Filter 4. remove sequences of bad length

lengths <- nchar(as.vector(filtered_data$seq))
filtered_data <- filtered_data[(lengths-313)%%3 == 0,]
lengths <- nchar(as.vector(filtered_data$seq))

print('Filter 4 finished')

#####
print('starting Filter 5: remove numts')
# Filter 5. remove numts

no_numts_data <- c()
numts_seqs <- c()

number_of_motus <- length(unique(filtered_data$motu))
for (i in 1:number_of_motus) {
  motu <- unique(filtered_data$motu)[i]
  datas <- filtered_data[filtered_data$motu==motu,]
  is_metazoa <- motu_taxa$Metazoa[motu_taxa$id==as.character(motu)]
  datas_length <- nchar(as.vector(datas$seq))
  newlist <- numts(datas, is_metazoa = is_metazoa, motu = motu, datas_length = datas_length)
  no_numts_data <- rbind(no_numts_data,newlist[['no_numts_data']])
  numts_seqs <- rbind(numts_seqs,newlist[['numts_seqs']])
  print(paste0((i/number_of_motus*100),'% run'))
}

filtered_data <- no_numts_data

print('Filter 5 finished')

print('writing outputs')

# create MOTU table with taxo
original_data <- original_data[match(unique(filtered_data$motu),original_data$id),]
motu_abund <- lapply(original_data$id, FUN = function(x){
  data_motu <- filtered_data[as.character(filtered_data$motu) == x ,grep('2018',colnames(filtered_data))]
  if (dim(data_motu)[1]>1) {
    return(colSums(data_motu))
  } else { 
    return(data_motu)
    }
} )

original_data <- cbind(original_data,bind_rows(lapply(motu_abund, as.data.frame.list)))


write.csv(filtered_data, paste(output_dir,"ESV_final_table.csv", sep = ""), row.names = F)
write.csv(numts_seqs, paste(output_dir,"ESV_numts.csv", sep = ""), row.names = F)
write.csv(original_data, paste(output_dir,"MOTU_final_table.csv", sep = ""), row.names = F)

print('finished')
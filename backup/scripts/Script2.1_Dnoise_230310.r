library("optparse")
library("stringr")
library("dplyr")
source('scripts/Script2.2_numts.R')

output_dir = "tmp/finalfiles"

print('loading data')

ASV_data_initial <- read.csv("tmp/ASV_Adcorr.csv")

# Load file with MJOLNIR agnomens and corresponding original names
metadata <- read.table("results/metadata/COSQ_metadata_new.tsv",sep="\t",head=T,stringsAsFactors = F)
#sample_metadata_sorted <- read.table(metadata2) # this would be for the metadata in an sorted way
original_data <- read.csv("results/COSQ_Curated_LULU.tsv", sep='\t')
# Are only metazoans used finally? Would like to include other groups as well
motu_taxa <- data.frame('id' = original_data$id, 'Metazoa' = c(original_data$kingdom == 'Metazoa' & !is.na(original_data$kingdom)))

print('data loaded')

#Select sample columns only
sample_cols <- grep("sample",colnames(ASV_data_initial))
#Remove the prefix "sample." from sample names
sample_names <- gsub("sample\\.","",names(ASV_data_initial[sample_cols]))

# Change agnomens to original names
new_sample_names <- metadata$original_samples[match(sample_names,metadata$mjolnir_agnomens)]
# The below did not reveal any NAs, so did not run
#new_sample_names[is.na(new_sample_names)] <- gsub("^","EMPTY",as.character(c(1:sum(is.na(new_sample_names_ASV)))))
names(ASV_data_initial)[sample_cols] <- new_sample_names

original_data <- original_data[,-grep('sample',colnames(original_data))]

neg_samples <- grep("CNE|SN|WN|control", colnames(ASV_data_initial))
neg_data <- ASV_data_initial[,neg_samples]

# Remove negative controls from data
ASV_data_initial <- ASV_data_initial[,-neg_samples]

# Merge data from the same sample but different sequencing runs
## First rename row names to ASV id
row.names(ASV_data_initial)<-ASV_data_initial$id

## Then, create OTU table
### First check where the first sample column is
ASV_data_initial[1,1:5]
### Check where the last sample column is
n<-ncol(ASV_data_initial)
ASV_data_initial[1,(n-1):n]
### Extract all sample columns
otu <- ASV_data_initial[,5:(n-1)] 
### Transpose
otu_t <- t(otu) 
otu_t <- data.frame(otu_t)
### Make a new column containing the sample name and PCR replicate number (without sequencing run)
otu_t$root <- factor(substr(row.names(otu_t),1,nchar(row.names(otu_t))-2))
### Aggregate data from the same sample and PCR replicate
otu_agg <- group_by(otu_t, root) %>% summarise(across(everything(), list(sum)))
otu_agg <- as.data.frame(otu_agg)
ncol(otu_agg)
### Remove suffix created by dplyr
otu_agg <- otu_agg %>% rename_with(~str_remove(., '_1'))
### Check that no columns were accidentally removed
ncol(otu_agg)
row.names(otu_agg)<-otu_agg$root
otu_agg<-within(otu_agg,rm("root"))
otu_agg_t <- t(otu_agg)
read_data <- data.frame(otu_agg_t,check.names=F)

# Merge new OTU table with taxonomy and sequences
initial_data <- ASV_data_initial[,c(1,2,4)]

seq_data <- ASV_data_initial$sequence

merged_data <- cbind(initial_data, read_data, seq_data)

#rm(ASV_data_initial, initial_data, seq_data)

#####
# print('starting Filter 1: remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance remove blank and NEG samples')
# # Filter 1. remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance
# # remove blank and NEG samples
# data_neg_filt_deleted <- merged_data[rowSums(neg_data)/rowSums(cbind(merged_data[,grep("2018",colnames(merged_data))],neg_data)) > 0.1,]
# merged_data <- merged_data[rowSums(neg_data)/rowSums(cbind(merged_data[,grep("2018",colnames(merged_data))],neg_data)) <= 0.1,]
# 
# # Add samples names and write data into neg_filtrates directory
# write.csv(merged_data, file = paste0(output_dir,"ASV_negfilt.csv"),row.names = F)
# 
# write.csv(data_neg_filt_deleted, file = paste0(output_dir,"ASV_negfilt_deleted.csv"),row.names = F)
# rm(data_neg_filt_deleted)
# 
# print('Filter 1 finished')

#####
# Filter 2. Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below X% of the total reads of this sample

min_relative <- 0.001
relabund <- function(x,min_relative) if (sum(x)>0) x/sum(x) < min_relative else FALSE

change_matrix <- do.call("cbind",apply(merged_data[,colnames(read_data)], 2, relabund, min_relative=min_relative)) & merged_data[,colnames(read_data)]>0

relabund_changed <- data.frame(ASV_id_modified = rownames(change_matrix[rowSums(change_matrix)>0,]),
                               samples = vapply(rownames(change_matrix[rowSums(change_matrix)>0,]), function(x,change_matrix){
                                 return(paste(colnames(change_matrix)[change_matrix[rownames(change_matrix)==x,]],collapse = "|"))
                               }, FUN.VALUE = "string", change_matrix = change_matrix))
merged_data[,colnames(read_data)][change_matrix] <- 0

merged_data$COUNT <- rowSums(merged_data[,colnames(read_data)])

min_reads <- 5

message("RAGNAROC is removing MOTUs with less than ",min_reads," total reads.")
merged_data <- merged_data[merged_data$COUNT >= min_reads,]

# remove numts
message("numts will be removed")
no_ASV_before_numts <- dim(merged_data)[1]
# First, remove sequences of incorrect length
lengths <- nchar(as.vector(merged_data$seq_data))
merged_data <- merged_data[(lengths-313)%%3 == 0,]
lengths <- nchar(as.vector(merged_data$seq_data))

no_numts_data <- c()
numts_seqs <- c()

number_of_motus <- length(unique(merged_data$motu))

cores <- 1
numts_ASV <- parallel::mclapply(1:number_of_motus,function(i,merged_data,motu_taxa){
  motu <- unique(merged_data$motu)[i]
  datas <- merged_data[merged_data$motu==motu,]
  is_metazoa <- motu_taxa$Metazoa[motu_taxa$id==as.character(motu)]
  datas_length <- nchar(as.vector(datas$seq_data))
  newlist <- numts(datas, is_metazoa = is_metazoa, motu = motu, datas_length = datas_length)
  return(newlist)
 },merged_data=merged_data,motu_taxa=motu_taxa,mc.cores = cores)
#numts_ASV_all <- do.call("rbind",numts_ASV)
ID_numts_ASV = unlist( lapply(numts_ASV, function(x) x$no_numts_data$id) )

final_data <- merged_data[merged_data$id %in% ID_numts_ASV,]

before<-length(merged_data$id)
after <- length(ID_numts_ASV)
no_numts <- before - after

message("",no_numts," numts removed. ",before," ASVs was reduced to ",after," ASVs")

# Remove old count column
merged_data = subset(merged_data, select = -c(count))

write.table(final_data,"results/COSQ_final_ASV.tsv",row.names = F,sep="\t",quote = F)
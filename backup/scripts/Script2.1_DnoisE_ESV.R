library("optparse")
library("stringr")
library("dplyr")
source('scripts/Script2.2_numts.R')

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="character"),
  make_option(c("-a", "--metadata"), type="character", default=NULL, metavar="character"),
  make_option(c("-b", "--original_data"), type="character", default=NULL, metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("input_file needed", call.=FALSE)
}
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("metadata needed", call.=FALSE)
}
if (is.null(opt$original_data)){
  print_help(opt_parser)
  stop("original_data needed", call.=FALSE)
}
if (is.null(opt$output_dir)){
  print_help(opt_parser)
  stop("output_dir needed", call.=FALSE)
}


# metadata <- read.csv('../../dades_originals/data_precuration/PHY1_sample_metadata.csv', stringsAsFactors = F)
# sample_metadata_sorted <- read.csv('../../Analisis/inputfiles/PHY_sample_metadata_ordered.csv', stringsAsFactors = F)

# ESV_data_initial <- read.csv('~/Documents/ESV_Adcorr.csv')
# load('metadata.RData')
# original_data <- read.csv("../../Pipeline/PHY1.Curated_LULU_euk.csv", sep='\t')
# motu_taxa <- data.frame('id' = original_data$id, 'Metazoa' = c(original_data$kingdom_name == 'Metazoa'))

output_dir = "tmp/finalfiles"
input = "tmp/ESV_Adcorr.csv"
metadata = "results/metadata/COSQ_metadata.tsv"

print('loading data')

#ESV_data_initial <- read.csv(opt$input)
ESV_data_initial <- read.csv(input)
load(metadata)
original_data <- read.csv(opt$original_data, sep='\t')
motu_taxa <- data.frame('id' = original_data$id, 'Metazoa' = c(original_data$kingdom_name == 'Metazoa'))

print('data loaded')

samples_col <- grep("sample",colnames(ESV_data_initial))
# colnames(ESV_data_initial)[samples_col] <- sample_metadata_NEG4B$codi
colnames(ESV_data_initial)[samples_col] <- sample_metadata$codi
original_data <- original_data[,-grep('sample',colnames(original_data))]

neg_samples <- grep("NEG|blank", colnames(ESV_data_initial))
neg_data <- ESV_data_initial[,neg_samples]

ESV_data_initial <- ESV_data_initial[,-neg_samples]
#samples_col <- grep("_2018",colnames(ESV_data_initial))

initial_data <- ESV_data_initial[,c(1,2,4)]

ordered_data <- ESV_data_initial[,samples_col] %>% select(sample_metadata_sorted$codi)

seq_data <- ESV_data_initial$sequence

filtered_data <- cbind(initial_data, ordered_data, seq_data)

rm(ESV_data_initial, initial_data, ordered_data, seq_data)


#####
print('starting Filter 1: remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance remove blank and NEG samples')
# Filter 1. remove any MOTU for which abundance in the blank or negative controls was higher than 10% of its total read abundance
# remove blank and NEG samples
data_neg_filt_deleted <- filtered_data[rowSums(neg_data)/rowSums(cbind(filtered_data[,grep("2018",colnames(filtered_data))],neg_data)) > 0.1,]
filtered_data <- filtered_data[rowSums(neg_data)/rowSums(cbind(filtered_data[,grep("2018",colnames(filtered_data))],neg_data)) <= 0.1,]

# Add samples names and write data into neg_filtrates directory
write.csv(filtered_data, file = paste0(output_dir,"ESV_negfilt.csv"),row.names = F)

write.csv(data_neg_filt_deleted, file = paste0(output_dir,"ESV_negfilt_deleted.csv"),row.names = F)
rm(data_neg_filt_deleted)

print('Filter 1 finished')


#####
print('starting Filter 2: Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below 0.005% of the total reads of this sample')
# Filter 2. Apply a minimum relative abundance threshold for each sample, setting to zero any abundance below 0.005% of the total reads of this sample

modified_ESV <- list()

relabund <- function(x){
  x/sum(x)
}



for (i in sample_metadata_sorted$codi) {
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

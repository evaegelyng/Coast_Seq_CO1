library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")

#Load tables

#This script should be run from the "results" folder, using the metabar_2021 environment. 

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("COSQ_otu_phyloseq.txt", sep="\t", header=T, row.names=1,check.names=F))

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

p_DADAwang = phyloseq(OTU)

rm(otu_mat)
gc()

#Load metadata

p_metadata<-read.table("metadata/COSQ_metadata_complete.txt", sep="\t", header=T)
p_metadata2<-p_metadata[,-5]
p_metadata2$sr_PCR_r<-paste(p_metadata2$
sample_root, p_metadata2$PCR_replicate, sep="_")
sampledata = sample_data(data.frame(p_metadata2, row.names=p_metadata2$sample_ID, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_DADAwang, sampledata)
DADAwang1 = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
DADAwang1

#rm(p_DADAwang) This should not be necessary - Marcelo agrees
gc()

##Merging samples from different runs (separately according to substrate - memory use!)
p_metadata3<-p_metadata2[,c("sr_PCR_r","substrate_type","sample_ID")]
vyt = sample_data(data.frame(p_metadata3, row.names=p_metadata3$sample_ID, stringsAsFactors=FALSE))

sample_data(DADAwang1)<-vyt
sample_data(DADAwang1)$dummy_var <- 1

sddd<-subset_samples(DADAwang1, substrate_type=="sediment")
wddd<-subset_samples(DADAwang1, substrate_type=="water")

final_table_sedi<-merge_samples(sddd, "sr_PCR_r", fun=sum)
final_table_wate<-merge_samples(wddd, "sr_PCR_r", fun=sum)

new_s_otu<-as.matrix(data.frame(otu_table(final_table_sedi), check.names=F))
new_w_otu<-as.matrix(data.frame(otu_table(final_table_wate), check.names=F))

new_f_matrix<-rbind(new_s_otu, new_w_otu)
nrow(new_f_matrix)==nrow(new_s_otu)+nrow(new_w_otu)
nrow(new_f_matrix)==length(unique(p_metadata2$sr_PCR_r))
rownames(new_f_matrix)==unique(p_metadata2$sr_PCR_r)

rm(sddd,wddd)
gc()


############Fix metadata --- (re-build it from 0)

sample_ID<-rownames(new_f_matrix)

###Get sample source
source<-NA
metadata<-data.frame(cbind(sample_ID, source))
rownames(metadata)<-metadata$sample_ID
metadata$p_source<- sapply(strsplit(as.character(metadata$sample_ID), "_"), head, 1)
metadata$ppn<-as.character(gsub('\\d','', metadata$p_source))

metadata$source<-ifelse(metadata$p_source=="CNE"|metadata$p_source=="2CCNE", 
as.character("CNE"), ifelse(metadata$p_source=="control", as.character("control"), ifelse(grepl("SN", metadata$ppn, fixed=T)|grepl("WN", metadata$ppn, fixed=T), as.character("NTC"), as.character("Field_sample"))))

metadata$root<-gsub("_+[^_]+$", "",metadata$sample_ID)
metadata$po<- sapply(strsplit(as.character(metadata$sample_ID), "_"), tail, 1)
metadata$PCR_replicate<-as.integer(gsub("\\D", "", metadata$po))
metadata$pn<-gsub('\\D','_', metadata$sample_ID)
metadata$pn2<-gsub(".*_(.+)__.*", "\\1", metadata$pn)
metadata$pn3<-gsub(".*[C]([^.]+)[_].*", "\\1", metadata$sample_ID)

metadata$cluster<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn2), ifelse(metadata$source=="control", as.integer(gsub('\\D','', metadata$pn3)), NA))

metadata$pn6<- sapply(strsplit(as.character(metadata$pn), "__"), tail, 1)
metadata$pn7<-sapply(strsplit(as.character(metadata$pn6), "_"), head, 1)
metadata$field_replicate<-ifelse(metadata$source=="Field_sample", as.integer(metadata$pn7), NA)
metadata$habitat<-ifelse(metadata$source=="Field_sample", as.character(gsub('\\d','', metadata$pn3)), NA)

extraction_refs<-read.table("~/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/backup/data/raw_data/extraction_refs_both_seasons.txt", sep="\t", header=T, row.names=1)
extraction_refs$sample_root<-row.names(extraction_refs)
extraction_refs$extraction_refs<-sapply(strsplit(as.character(extraction_refs$extraction_number), "_"), head, 1)

metadata$extraction_refs<-extraction_refs$extraction_refs[match(metadata$root, extraction_refs$sample_root)]

PSU_refs<-read.table("metadata/COSQ_metadata_reps.tsv", sep="\t", header=T)
metadata$PSU_refs<-PSU_refs$PSU[match(metadata$root, PSU_refs$root)]
metadata$substrate_type<-PSU_refs$substrate_type[match(metadata$root, PSU_refs$root)]
metadata$season<-PSU_refs$season[match(metadata$root, PSU_refs$root)]

new_metadata<-metadata[,c("sample_ID", "root", "source", "season", "substrate_type", "cluster", "habitat", "field_replicate", "extraction_refs", "PSU_refs", "PCR_replicate")]


###Setup new phyloseq object
new_sampledata = sample_data(data.frame(new_metadata, row.names=new_metadata$sample_ID, stringsAsFactors=FALSE))

t_new_f_matrix<-t(new_f_matrix)
OTU2 = otu_table(t_new_f_matrix, taxa_are_rows = TRUE)

final_table = phyloseq(OTU2)

final_table2 = merge_phyloseq(final_table, new_sampledata)

# ? otu_table(final_table)<-t(otu_table(final_table))
sum(sample_sums(DADAwang1))==sum(sample_sums(final_table2))

rm(DADAwang1)

##backup main table
str(t_new_f_matrix)

new_otu_mat<-as.matrix(data.frame(otu_table(t_new_f_matrix, taxa_are_rows = TRUE), check.names=F))
str(new_otu_mat)

utus<-new_otu_mat
sum(colSums(utus))

###SPRING

dir.create("cleanup_ASV_wise")
dir.create("cleanup_ASV_wise/spring")

final_table_s<-subset_samples(final_table2, season=="spring")

#####

#Sed decon by PSU set
cat("spring_sed_NTC")
#Subset substrate type (sediment)
data_sedi<-subset_samples(final_table_s, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate. EES: Added commands to skip errors due to 
#some combinations of PSU and PCR replicates containing no NTCs (all sequences removed by MJOLNIR pipeline)  
  skip_to_next <- FALSE
  tryCatch(subset_samples(data_sedi, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr]), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }
  contam<-subset_samples(data_sedi, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)  
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs 
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/spring/cont_list_sed_ntc_spring",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(new_otu_mat))-sum(colSums(utus))
after_ntc_decon_sed_spring<-utus

tmp<-"../tmp/"
mandss<-as.matrix(after_ntc_decon_sed_spring)
saveRDS(mandss, paste(tmp,"mandss.rds",sep=""))

rm(new_otu_mat,mandss)
gc()

#Water decon by PSU set
cat("spring_wat_NTC")
#Subset substrate type (water)
data_wate<-subset_samples(final_table_s, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate. EES: Added commands to skip errors due to 
#some combinations of PSU and PCR replicates containing no NTCs (all sequences removed by MJOLNIR pipeline)  
  skip_to_next <- FALSE
  tryCatch(subset_samples(data_wate, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr]), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }  
  contam<-subset_samples(data_wate, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/spring/cont_list_wat_ntc_spring",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_ntc_decon_sed_spring))-sum(colSums(utus))
after_ntc_decon_wat_spring<-utus

mandws<-as.matrix(after_ntc_decon_wat_spring)
saveRDS(mandws, paste(tmp,"mandws.rds",sep=""))

rm(mandws)
gc()

#Sed decon by extraction set
cat("spring_sed_CNE")
#Subset substrate type (sediment)
otu_table(final_table_s) = otu_table(after_ntc_decon_sed_spring, taxa_are_rows = TRUE)
pdata_sedi<-subset_samples(final_table_s, substrate_type=="sediment")
#Check whether any extraction sets are missing a CNE
samples<-pdata_sedi@sam_data[which(pdata_sedi@sam_data$source=="Field_sample"),]
levels(factor(samples$extraction_refs))
cnes<-pdata_sedi@sam_data[which(pdata_sedi@sam_data$source=="CNE"),]
levels(factor(cnes$extraction_refs))
#Exclude extraction sets with no CNEs in the final dataset
data_sedi<-subset_samples(pdata_sedi, !extraction_refs=="S2" & !extraction_refs=="S4")

tm<-data.frame(sample_data(data_sedi))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_sedi, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/spring/cont_list_sed_extr_spring",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ntc_decon_wat_spring))-sum(colSums(utus))
after_ext_decon_sed_spring<-utus

maedss<-as.matrix(after_ext_decon_sed_spring)
saveRDS(maedss, paste(tmp,"maedss.rds",sep=""))

rm(after_ntc_decon_sed_spring,maedss)
gc()

#Wat decon by extraction set
cat("spring_wat_CNE")
#Subset substrate type (water)
otu_table(final_table_s) = otu_table(after_ntc_decon_wat_spring, taxa_are_rows = TRUE)
pdata_wate<-subset_samples(final_table_s, substrate_type=="water")
#Check whether any extraction sets are missing a CNE
samples<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="Field_sample"),]
levels(factor(samples$extraction_refs))
cnes<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="CNE"),]
levels(factor(cnes$extraction_refs))
#Exclude extraction sets with no CNEs in the final dataset
#data_wate<-subset_samples(pdata_wate, !extraction_refs=="XX") # No CNEs were missing
data_wate<-pdata_wate

tm<-data.frame(sample_data(data_wate))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_wate, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/spring/cont_list_wat_extr_spring",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_sed_spring))-sum(colSums(utus))
after_ext_decon_wate_spring<-utus

maedws<-as.matrix(after_ext_decon_wate_spring)
saveRDS(maedws, paste(tmp,"maedws.rds",sep=""))

rm(after_ntc_decon_wat_spring,after_ext_decon_sed_spring,maedws)
gc()

#Wat decon by field control
cat("spring_wat_field_control")
#Subset substrate type (water)
otu_table(final_table_s) = otu_table(after_ext_decon_wate_spring, taxa_are_rows = TRUE)
pdata_wate<-subset_samples(final_table_s, substrate_type=="water")
####exclude clusters with no field control in the final dataset
controls<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="control"),]
list<-levels(factor(controls$cluster))
idx= pdata_wate@sam_data$cluster %in% as.numeric(list)
data_wate<-subset_samples(pdata_wate, idx)
tm<-data.frame(sample_data(data_wate))
clst<-levels(as.factor(tm[,"cluster"]))

for (ps in 1:length(clst)){
#summarize contaminants per cluster  
  contam<-subset_samples(data_wate, source=="control"&cluster==clst[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(clst[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2))
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, source=="Field_sample"&cluster==clst[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2))
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(clst[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/spring/cont_list_wat_clst_spring",clst[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_wate_spring))-sum(colSums(utus))
after_clst_decon_wate_spring<-utus

macdws<-as.matrix(after_clst_decon_wate_spring)
saveRDS(macdws, paste(tmp,"macdws.rds",sep=""))

rm(after_ext_decon_wate_spring,macdws)
gc()


#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#>>>>>>>>><<<<<<<<<>>>>>>>>>><<<<<<<<#
#Season2 - keep fixing same utus object#

###AUTUMN

dir.create("cleanup_ASV_wise/autumn")

final_table_a<-subset_samples(final_table2, season=="autumn")

#Sed decon by PSU set
cat("autumn_sed_NTC")
#Subset substrate type (sediment)
data_sedi<-subset_samples(final_table_a, substrate_type=="sediment")
tm<-data.frame(sample_data(data_sedi))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate. EES: Added commands to skip errors due to 
#some combinations of PSU and PCR replicates containing no NTCs (all sequences removed by MJOLNIR pipeline)  
  skip_to_next <- FALSE
  tryCatch(subset_samples(data_sedi, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr]), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }  
  contam<-subset_samples(data_sedi, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs 
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/autumn/cont_list_sed_ntc_autumn",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_clst_decon_wate_spring))-sum(colSums(utus))
after_ntc_decon_sed_autumn<-utus

mandsa<-as.matrix(after_ntc_decon_sed_autumn)
saveRDS(mandsa, paste(tmp,"mandsa.rds",sep=""))

rm(after_clst_decon_wate_spring,mandsa)
gc()

#Water decon by PSU set
cat("autumn_wat_NTC")
#Subset substrate type (water)
data_wate<-subset_samples(final_table_a, substrate_type=="water")
tm<-data.frame(sample_data(data_wate))
psus<-levels(as.factor(tm[,"PSU_refs"]))
prep<-levels(as.factor(tm[,"PCR_replicate"]))

for (ps in 1:length(psus)){
for (pr in 1:length(prep)){
#summarize contaminants per PSU&PCRreplicate. EES: Added commands to skip errors due to 
#some combinations of PSU and PCR replicates containing no NTCs (all sequences removed by MJOLNIR pipeline)  
  skip_to_next <- FALSE
  tryCatch(subset_samples(data_wate, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr]), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) { next }  
  contam<-subset_samples(data_wate, source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  if(sum(colSums(otu_table(contam)))==0) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="NTC"&PSU_refs==psus[ps]&PCR_replicate==prep[pr])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")],   ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(pr)
  cat(psus[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/autumn/cont_list_wat_ntc_autumn",psus[ps],prep[pr],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }

}
}

sum(colSums(after_ntc_decon_sed_autumn))-sum(colSums(utus))
after_ntc_decon_wat_autumn<-utus

mandwa<-as.matrix(after_ntc_decon_wat_autumn)
saveRDS(mandwa, paste(tmp,"mandwa.rds",sep=""))

rm(mandwa)
gc()

#Sed decon by extraction set
cat("autumn_sed_CNE")
#Subset substrate type (sediment)
otu_table(final_table_a) = otu_table(after_ntc_decon_sed_autumn, taxa_are_rows = TRUE)
pdata_sedi<-subset_samples(final_table_a, substrate_type=="sediment")
#Check whether any extraction sets are missing a CNE
samples<-pdata_sedi@sam_data[which(pdata_sedi@sam_data$source=="Field_sample"),]
levels(factor(samples$extraction_refs))
cnes<-pdata_sedi@sam_data[which(pdata_sedi@sam_data$source=="CNE"),]
levels(factor(cnes$extraction_refs))
#Exclude extraction sets with no CNEs in the final dataset
#data_sedi<-subset_samples(pdata_sedi, !extraction_refs=="2CVS4") # No CNEs were missing
data_sedi<-pdata_sedi

tm<-data.frame(sample_data(data_sedi))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction  
  contam<-subset_samples(data_sedi, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_sedi, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/autumn/cont_list_sed_extr_autumn",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ntc_decon_wat_autumn))-sum(colSums(utus))
after_ext_decon_sed_autumn<-utus

maedsa<-as.matrix(after_ext_decon_sed_autumn)
saveRDS(maedsa, paste(tmp,"maedsa.rds",sep=""))

rm(after_ntc_decon_sed_autumn,maedsa)
gc()

#Wat decon by extraction set
cat("autumn_wat_CNE")
#Subset substrate type (water)
otu_table(final_table_a) = otu_table(after_ntc_decon_wat_autumn, taxa_are_rows = TRUE)
pdata_wate<-subset_samples(final_table_a, substrate_type=="water")
#Check whether any extraction sets are missing a CNE
samples<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="Field_sample"),]
levels(factor(samples$extraction_refs))
cnes<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="CNE"),]
levels(factor(cnes$extraction_refs))
#Exclude extraction sets with no CNEs in the final dataset
data_wate<-subset_samples(pdata_wate, !extraction_refs=="MPAA15")

tm<-data.frame(sample_data(data_wate))
extr<-levels(as.factor(tm[,"extraction_refs"]))

for (ps in 1:length(extr)){
#summarize contaminants per extraction.
  contam<-subset_samples(data_wate, source=="CNE"&extraction_refs==extr[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(extr[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, !source=="CNE"&extraction_refs==extr[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(extr[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/autumn/cont_list_wat_extr_autumn",extr[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_sed_autumn))-sum(colSums(utus))
after_ext_decon_wate_autumn<-utus

maedwa<-as.matrix(after_ext_decon_wate_autumn)
saveRDS(maedwa, paste(tmp,"maedwa.rds",sep=""))

rm(after_ext_decon_sed_autumn,after_ntc_decon_wat_autumn,maedwa)
gc()

#Wat decon by field control
cat("autumn_wat_field_control")
#Subset substrate type (water)
otu_table(final_table_a) = otu_table(after_ext_decon_wate_autumn, taxa_are_rows = TRUE)
pdata_wate<-subset_samples(final_table_a, substrate_type=="water")
####exclude clusters with no field control in the final dataset
controls<-pdata_wate@sam_data[which(pdata_wate@sam_data$source=="control"),]
list<-levels(factor(controls$cluster))
idx= pdata_wate@sam_data$cluster %in% as.numeric(list)
data_wate<-subset_samples(pdata_wate, idx)


tm<-data.frame(sample_data(data_wate))
clst<-levels(as.factor(tm[,"cluster"]))

for (ps in 1:length(clst)){
#summarize contaminants per cluster  
  contam<-subset_samples(data_wate, source=="control"&cluster==clst[ps])  
  if(sum(colSums(otu_table(contam)))==0) next
  cat(clst[ps])
  cat("\n")
################################################################# 
  p_contam2 = filter_taxa(contam, function(x) sum(x) > 0, TRUE)
  contam2 = prune_samples(sample_sums(p_contam2)>0,p_contam2)
  ag_contam2<-data.frame(otu_table(contam2),check.names=F)
  ag_contam2[, "max_control"] <- apply(ag_contam2, 1, max)
  ag_contam2[, "sum_control"] <- apply(ag_contam2, 1, function(x) sum(x)-max(x))
  contam_list<-taxa_names(contam2)
#summarize contaminants in field samples  
  ncontam<-subset_samples(data_wate, source=="Field_sample"&cluster==clst[ps])
  ncontam2 = prune_taxa(contam_list, ncontam)
  ag_ncontam2<-data.frame(otu_table(ncontam2),check.names=F)
  r_ag_ncontam2<-ag_ncontam2
  ag_ncontam2[, "max_field"] <- apply(ag_ncontam2, 1, max)
  ag_ncontam2[, "sum_field"] <- apply(ag_ncontam2, 1, function(x) sum(x)-max(x))
#merge the controls and field samples' summaries  
  s_ct<-cbind(taxa_names(contam2), ag_ncontam2[, c("max_field", "sum_field")], ag_contam2[, c("max_control", "sum_control")])
###################################################
  if(length(unique(s_ct$max_control-s_ct$max_field>0))==1 && unique(s_ct$max_control-s_ct$max_field>0)==FALSE) next
  cat(clst[ps])
  cat("\n")
  er_list<-rownames(subset(s_ct,s_ct$max_control-s_ct$max_field>0))
#save list of contaminant OTUs  
  write.table(s_ct[,-1], file=paste("cleanup_ASV_wise/autumn/cont_list_wat_clst_autumn",clst[ps],".txt", sep="_"), quote=FALSE, sep='\t', row.names=TRUE)
#Fix otu_table  (set 0 for OTUs with control counts > field samples counts)
  for (e in 1:length(er_list))
  {
    for (i in 1:ncol(r_ag_ncontam2))
    {
      utus[er_list[e], colnames(r_ag_ncontam2)[i]]<-0  
    }
  }
}

sum(colSums(after_ext_decon_wate_autumn))-sum(colSums(utus))
after_clst_decon_wate_autumn<-utus

macdwa<-as.matrix(after_clst_decon_wate_autumn)
saveRDS(macdwa, paste(tmp,"macdwa.rds",sep=""))

rm(after_clst_decon_wate_autumn,after_ext_decon_wate_autumn,macdwa)
rm(final_table_a,final_table_s,data_wate,data_sedi)
gc()

###Write cleaned otu table

write.table(utus, "cleaned_otu_table_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)

#
rm(utus)
gc()

###Write cleaned metadata

write.table(data.frame(sample_data(final_table2), check.names=F), "metadata/merged_seq_run_cleaned_metadata_ASV_wise.txt", sep="\t", quote=FALSE, row.names=TRUE)

###Summarize and write new ASV table

dir.create("cleanup_ASV_wise/summary")
summary<-"cleanup_ASV_wise/summary/"

new_otu_ss<-subset_samples(final_table2, season=="spring")
new_otu_aa<-subset_samples(final_table2, season=="autumn")

no<-otu_table(final_table2)

#
rm(final_table2)
gc()

nd<-data.frame(no, check.names=F)
new_m_otu<-as.matrix(nd)

nos<-otu_table(new_otu_ss)
nds<-data.frame(nos, check.names=F)
new_m_otu_s<-as.matrix(nds)

noa<-otu_table(new_otu_aa)
nda<-data.frame(noa, check.names=F)
new_m_otu_a<-as.matrix(nda)

rm(no, nd)
gc()

rm(nos, nds)
gc()

rm(noa, nda)
gc()


mandss<-readRDS("../tmp/mandss.rds")
mandws<-readRDS("../tmp/mandws.rds")
maedss<-readRDS("../tmp/maedss.rds")
maedws<-readRDS("../tmp/maedws.rds")
macdws<-readRDS("../tmp/macdws.rds")
mandsa<-readRDS("../tmp/mandsa.rds")
mandwa<-readRDS("../tmp/mandwa.rds")
maedsa<-readRDS("../tmp/maedsa.rds")
maedwa<-readRDS("../tmp/maedwa.rds")
macdwa<-readRDS("../tmp/macdwa.rds")

andss<-mandss[,colnames(new_m_otu_s)]
andws<-mandws[,colnames(new_m_otu_s)]
aedss<-maedss[,colnames(new_m_otu_s)]
aedws<-maedws[,colnames(new_m_otu_s)]
acdws<-macdws[,colnames(new_m_otu_s)]

andsa<-mandsa[,colnames(new_m_otu_a)]
andwa<-mandwa[,colnames(new_m_otu_a)]
aedsa<-maedsa[,colnames(new_m_otu_a)]
aedwa<-maedwa[,colnames(new_m_otu_a)]
acdwa<-macdwa[,colnames(new_m_otu_a)]

scus1<-colSums(new_m_otu_s)
scus2<-colSums(andss)
scus3<-colSums(andws)
scus4<-colSums(aedss)
scus5<-colSums(aedws)
scus6<-colSums(acdws)

summary_clean_up_spring1<-cbind(scus1,scus2)
summary_clean_up_spring2<-cbind(scus3,scus4)
summary_clean_up_spring3<-cbind(scus5,scus6)

rm(new_m_otu_s)
gc()

summary_clean_up_spring<-cbind(summary_clean_up_spring1,summary_clean_up_spring2,summary_clean_up_spring3)
colnames(summary_clean_up_spring)<-c("raw","ntc_sed","ntc_wat","cne_sed","cne_wat","cntrl_wat")
write.table(summary_clean_up_spring,paste(summary,"summary_clean_up_spring.txt",sep=""))

scua1<-colSums(new_m_otu_a)
scua2<-colSums(andsa)
scua3<-colSums(andwa)
scua4<-colSums(aedsa)
scua5<-colSums(aedwa)
scua6<-colSums(acdwa)

summary_clean_up_autumn1<-cbind(scua1,scua2)
summary_clean_up_autumn2<-cbind(scua3,scua4)
summary_clean_up_autumn3<-cbind(scua5,scua6)

rm(new_m_otu_a)
gc()

summary_clean_up_autumn<-cbind(summary_clean_up_autumn1,summary_clean_up_autumn2,summary_clean_up_autumn3)
colnames(summary_clean_up_autumn)<-c("raw","ntc_sed","ntc_wat","cne_sed","cne_wat","cntrl_wat")
write.table(summary_clean_up_autumn,paste(summary,"summary_clean_up_autumn.txt",sep=""))


csn<-colSums(new_m_otu)
csn1<-colSums(mandss)
csn2<-colSums(mandws)
csn3<-colSums(maedss)
csn4<-colSums(maedws)
csn5<-colSums(macdws)
csn6<-colSums(mandsa)
csn7<-colSums(mandwa)
csn8<-colSums(maedsa)
csn9<-colSums(maedwa)
csn10<-colSums(macdwa)

scb<-cbind(csn,csn1)
scb1<-cbind(csn2,csn3)
scb2<-cbind(csn4,csn5)
scb3<-cbind(csn6,csn7)
scb4<-cbind(csn8,csn9)

scb_1<-cbind(scb,scb1)
scb1_1<-cbind(scb2,scb3)
scb2_1<-cbind(scb4,csn10)

rm(new_m_otu)
gc()

psummary_clean_up_both_seasons<-cbind(scb_1,scb1_1)
summary_clean_up_both_seasons<-cbind(psummary_clean_up_both_seasons,scb2_1)

colnames(summary_clean_up_both_seasons)<-c("raw","ntc_sed_s","ntc_wat_s","cne_sed_s","cne_wat_s","cntrl_wat_s","ntc_sed_a","ntc_wat_a","cne_sed_a","cne_wat_a","cntrl_wat_a")
write.table(summary_clean_up_both_seasons,paste(summary,"summary_clean_up_both_seasons.txt",sep=""))

###Test cleaned table

otu_matt<-as.matrix(read.table("cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

OTUt = otu_table(otu_matt, taxa_are_rows = TRUE)
p_DADAwangt = phyloseq(OTUt)

p_DADAwang
cat("total_reads_before_cleanup")
sum(sample_sums(p_DADAwang))
cat("\n")
p_DADAwangt
cat("total_reads_after_cleanup")
sum(sample_sums(p_DADAwangt))

## Count number of ASVs after cleaning
final<-read.table("COSQ_final_dataset.tsv", sep="\t", header=T, check.names=F)
otu_mat<-as.data.frame(otu_matt)
otu_mat$asvs<-final$cluster_weight[match(row.names(otu_mat),final$id)]
sum(otu_mat$asvs)
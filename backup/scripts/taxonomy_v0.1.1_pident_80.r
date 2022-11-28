#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
######################################################################################################################################################

# Taxonomical classification of ASVs

######################################################################################################################################################
#
# Authors: This script was mainly written by Tobias G. Frøslev (Copenhagen University). 
# It has here been modified by Eva Egelyng Sigsgaard (EES) such that instead of a fixed threshold ("upper margin") determining which BLAST hits to include for classification, this threshold
# is determined for each query sequence as the minimum similarity obtained for the best matching taxid. Adrián Gómez Repollés (AGR) and Caitlin Kim Frankish (CKF) have written the code to correct outdated taxids that are no longer valid, 
# as they have been merged with other taxids. The functions for indicating possible misidentifications of reference sequences were mainly developed by Emil Ellegaard Thomassen (EET).

# The script requires a BLAST output from a set of ASVs
# Minimum fields required are qseqid, staxid, pident, ssciname and evalue

# Instructions
# 1) Load the three functions below: assign_taxonomy, prefilter, get_classification, evaluate_classification
# 2) Classify your ASVs by running the wrapper function assign_taxonomy like this:
#   classified_table <- assign_taxonomy(INPUT.blasthits, lower_margin = 2, remove = c("unwanted_taxon1","unwanted_taxon2"))
#
# Explanation to input
#   INPUT.blasthits are the BLAST results
#   upper_margin is the margin used for suboptimal hits used for classification - e.g. a margin of 0.5 means that hits of 100% to 99.5% is used or 95% to 94.5%, if the best hit is 100% or 95%, respectively.
#   lower_margin: hits down to this margin from the best hit are shown in the output as alternative possibilities, but are not used for taxonomic classification.
#   remove: a vector of taxa to exclude from the evaluation. Could be e.g. remove = c("uncultured","environmental") to exclude hits with no precise annotation, or names of species known not to be in the study area.
#
# Explanation to output
#   The output is a list with
#   $classified_table : the table with all ASVs classified. One row per ASV
#        This table contains the estimated best classification at all taxonomic levels, based on the hits in the upper_margin (hits are weighted by evalue), 
#        Each taxonomic level gets a score indicating the agreement on the selected classification at that level.
#        Also a string of alternatives and their matches (%). This string includes hits from both upper and lower margin
#   $all_classifications: this is the table used to make the classified_table. It contains all hits above lower_margin for all ASVs and their classifications (only upper_margin).
#   ...and the input parameters

# Load required packages
library(taxizedb) # For retrieving taxonomic classifications
library(dplyr)
library(tidyr)
#library(stringr) # For splitting character strings into words. Only necessary when using BOLD+nt database (see below) 

# Provide API key for NCBI
options(ENTREZ_KEY="d48b9acccd829282ff4386f716341ed2f608")

# Read the completed BLAST results into a table
IDtable <- read.csv(file = args[1], sep='\t', header=F, as.is=TRUE)

# If using the GenBank nt database, use the following to add an empty column for ssciname. This is a temporary fix, as we have had problems retrieving scientific names from BLAST against our local reference database.
# If you would like to have the scientific names of each BLAST hit, try adding "ssciname" to the BLAST command in the gwf workflow file, and do not create the empty column.
IDtable$V16<-"NA"

# If using a combined BOLD+nt database built with the MARES pipeline, use the following commented lines to add sscinames
seqtab <- data.frame(readLines("/faststorage/project/eDNA/blastdb/Eukaryota_COI_NOBAR/taxid_process/Eukaryota_informative_name_table.tsv"))
seqtab$NWord <- sapply(strsplit(as.character(seqtab[,1]), " "), length)
colnames(seqtab)<-c("Name","NWord")

#Extract taxonomy from blast file and put in separate column called "taxa"
word<-stringr::word #Specifies that the function word should be taken from the stringr library
seqtab$Taxa <- ifelse(seqtab[,2] == 4,word(seqtab[,1],2,sep=" "),word(seqtab[,1],2,3,sep=" "))

#Extract TaxID from blast file and put in separate column called "TaxID"   
seqtab$TaxID <- ifelse(seqtab[,2] == 4,word(seqtab[,1],3,sep=" "),word(seqtab[,1],4,sep=" "))

#Add taxa column to IDtable based on matched TaxIDs
IDtable$V16 <- seqtab$Taxa[match(IDtable$V15,seqtab$TaxID)]
IDtable$V16<-gsub(" "," ",IDtable$V16) #This may be unnecessary

# Read the possible problematic taxids as a table
#MergedTaxIDs<-read.table("~/eDNA/faststorage/blastdb/nt_220203/taxdump/MergedTaxIDs", header=TRUE)
MergedTaxIDs<-read.table("~/eDNA/faststorage/blastdb/Euk_COI_NOBAR_2022/scripts/taxdump_edit/MergedTaxIDs", header=TRUE)

# Add header information
names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","staxid","ssciname")

# Extract only those rows where the query coverage is 100. First check whether there are in fact any hits with 100% query coverage.
    if (max(IDtable$qcovs) == 100 ) {
      IDtable <- IDtable[IDtable$qcovs==100,]  
    } else {
      readr::write_file("", args[2])
      readr::write_file("", args[3])  
      stop("Query coverage is less than 100% for all hits", call.=FALSE)
    }

# Extract only those rows where the pident is at least 90%. First check if the data contains hits with at least 90% sequence identity
#    if (max(IDtable$pident) >= 90 ) {
#      IDtable <- IDtable[IDtable$pident>=90,]  
#    } else {
#      readr::write_file("", args[2])
#      readr::write_file("", args[3])
#      stop("Sequence identity is less than 90% for all hits", call.=FALSE)
#    }

# Remove hits that are (definitely) not identified to species level (their scientific name only contains one word)
IDtable <- IDtable %>% mutate(species_level=if_else(sapply(strsplit(as.character(IDtable$ssciname)," "), length)==1, "no","yes"))
IDtable <- IDtable[IDtable$species_level=="yes",]

# The following is to define for each query sequence a minimum threshold of sequence similarity, which will determine whether a BLAST hit will be taken into account in the taxonomic classification of the query
# In the summary object below, the values needed to determine this threshold are calculated

# First, determine maximum and minimum similarity (pident) for each qseqid+taxid combination. Scientific names are unique for each taxid, and can thus also be included in the table.
summary<-do.call(data.frame,aggregate(pident~qseqid+staxid+ssciname,data=IDtable,FUN = function(x) c(max = max(x), min = min(x),n = length(x))))

# Sort the taxid hits by descending maximum similarity for each qseqid
summary<-summary[with(summary,order(summary$qseqid,-summary$pident.max)),]

# Add column indicating the maximum similarity for the best matching taxid(s)
summary$pident.max.best<-"NA" # Creating new column
for (i in unique (summary$qseqid)){
   summary[summary$qseqid==i,]$pident.max.best<-max(summary[summary$qseqid==i,]$pident.max)}   

# Add column for selecting unique qseqid+taxid combinations
summary$qseqid_staxid<-paste(summary$qseqid,summary$staxid,sep="_")

# Add a column indicating the best matching taxid(s)
summary$pident.best<-"NA" # Creating new column
for (i in unique (summary$qseqid_staxid)){
    if (summary[summary$qseqid_staxid==i,]$pident.max == summary[summary$qseqid_staxid==i,]$pident.max.best) {
      summary[summary$qseqid_staxid==i,]$pident.best<-"yes"
    } else {
      summary[summary$qseqid_staxid==i,]$pident.best<-"no"
    }
}

# For each qseqid, determine the minimum similarity for the best matching taxid(s) 
summary$pident.min.best<-"NA" # Creating new column
for (i in unique (summary$qseqid)){
   best<- summary[which(summary$qseqid==i & summary$pident.best=="yes"), ]
   pident.min.best<-min(best$pident.min)
   summary[summary$qseqid==i,]$pident.min.best<-pident.min.best}   

# Add a column that shows whether the present taxid overlaps in sequence similarity with the best matching taxid, and should therefore be included in the taxonomic classification
summary$include<-ifelse(summary$pident.max>=as.numeric(summary$pident.min.best),1,0)

# Add a column that indicates whether there is an unexpectedly large range of variation in sequence similarity within the same included taxid, suggesting possible misidentifications of specimens in the database. Here, the threshold is set at 3% variation, but the appropriate level depends on the metabarcode and target organisms.
summary$possible.misid.highrange<-"NA"  
for (i in unique (summary$qseqid_staxid)) {
  summary[summary$qseqid_staxid==i,]$possible.misid.highrange<-ifelse(summary[summary$qseqid_staxid==i,]$include==1 & (summary[summary$qseqid_staxid==i,]$pident.max - summary[summary$qseqid_staxid==i,]$pident.min) > 3,1,0)
}

# Test if one taxid in those to include stands out - adds 1 to possible.misid if pident.n of a certain taxid is 1 and the sum of all included sequences is greater than 2 times the number of taxids to include (suggests an underrepresented taxid) (Emil: 07/01/2022)
summary$possible.misid.outlier<-"NA"  
for (j in unique (summary$qseqid)) {
  for (i in unique (summary[summary$qseqid==j,]$qseqid_staxid)) {
    summary[summary$qseqid_staxid==i,]$possible.misid.outlier<-ifelse(summary[summary$qseqid_staxid==i,]$include==1 & summary[summary$qseqid_staxid==i,]$pident.n==1 & sum(summary[summary$qseqid==j & summary$include==1,]$pident.n)> 2*sum(summary[summary$qseqid==j,]$include==1),1,0)
  }
}

# Test if the identification is solely based on less than 3 sequences - to detect identifications because of 1 or 2 high-similarity hits, which could be erroneous (EET: 07/01/2022)
summary$possible.misid.few<-"NA"  
for (j in unique (summary$qseqid)) {
  for (i in unique (summary[summary$qseqid==j,]$qseqid_staxid)) {
    summary[summary$qseqid_staxid==i,]$possible.misid.few<-ifelse(summary[summary$qseqid_staxid==i,]$include==1 & sum(summary[summary$qseqid==j,]$include==1)==1 & summary[summary$qseqid_staxid==i,]$pident.n<3,1,0)
  }
}

#Add summary column to detect possible.misid based on at least one of the preceding tests (EET: 10/01/2022)
summary$possible.misid<-"NA"  
for (i in unique (summary$qseqid_staxid)) {
  summary[summary$qseqid_staxid==i,]$possible.misid<-ifelse(sum(as.integer(summary[summary$qseqid_staxid==i,]$possible.misid.highrange),as.integer(summary[summary$qseqid_staxid==i,]$possible.misid.outlier), as.integer(summary[summary$qseqid_staxid==i,]$possible.misid.few))>0,1,0)
}

# Write all the calculated values to a file
write.table(summary,file=args[2],sep="\t",row.names=FALSE)

# Add minimum similarity for the best matching taxid to IDtable from summary table
IDtable$pident.min.best<-"NA"
for (i in unique (IDtable$qseqid)){
   IDtable[IDtable$qseqid==i,]$pident.min.best<-summary[summary$qseqid==i,]$pident.min.best[1]
}

# Define the four taxonomy script functions (where FunctionX is a wrapper that runs 1, 2, 3 in one go):

# FunctionX
# Wrapper function using the three main functions - each step can also be done manually
assign_taxonomy <- function(table,lower_margin=2, remove = c("")) {       # EES removed constant upper margin specification
  pf <- prefilter(table, lower_margin, remove)   # EES removed constant upper margin specification
  
  # Replace old tax ids with new ones if they have a match in mergedtaxids dataframe #CKF 
  pf$OldTaxID<-as.numeric(pf$staxid) # Making an extra column on pf, same as staxid but with a new name so it can match with the mergedtaxid dataframe. EES added as.numeric to avoid error w. incompatible column classes
  pf_intermediate<-pf %>%
  dplyr::left_join(MergedTaxIDs, by=c("OldTaxID"))  # now we can join with MergedTaxIDs using the oldtaxid column
  
  # Print which IDs are being changed 
  pf_changed<-pf_intermediate %>%
  dplyr::filter(!is.na(NewTaxID))
  
  if (nrow(pf_changed)>0) {
  
  for (i in 1:nrow(pf_changed)) {
  pf_sub<-pf_changed[i,]
  print(paste("Warning:", pf_sub$qseqid, "has an outdated taxid. Overwriting outdated taxid:", pf_sub$OldTaxID, "to new taxid:", pf_sub$NewTaxID, sep=" "))
  }
  
 }  
  
  pf_new<-pf_intermediate %>%
  dplyr::mutate(staxid=ifelse(is.na(NewTaxID), staxid, NewTaxID)) %>% # If there is no match, then staxid stays the same, otherwise it takes the id in NewTaxID
  dplyr::select(-c(OldTaxID, NewTaxID)) # Now remove the two columns added to the dataframe, no longer needed 
  
  
  ##
  
  gc <- get_classification(pf_new)
  cf <- evaluate_classification(gc[[1]])     # AGR # gc[2] is a list of the taxids not classified 
  result <- list(classified_table=cf$taxonon_table, all_classifications=cf$all_taxa_table, all_classifications_summed=cf$all_taxa_table_summed, lower=lower_margin, removed=remove)  # EES removed upper_margin
  if (length(gc[[2]]) != 0) {
    print(paste0("Taxids not found in the classification: ", gc[2])) # AGR - Print the list of not matched taxids  
  }
  return(result) 
}

# Function1
# Filter data ASV-wise according to upper and lower margin set, and taxa to exclude
prefilter <- function(IDtable, lower_margin=2, remove = c("uncultured", "environmental")) {   # EES removed constant upper_margin specification
  new_IDtable <- IDtable[0,] # prepare filtered matchlist
  IDtable <- IDtable[!IDtable$staxid == "N/A",]
  ids <- names(table(IDtable$qseqid))
  i=1
  o=length(ids)
  for (name in ids) {
    test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
    if (nchar(remove[1])>0) {
      test2 <- test
      for (rm in 1:length(remove)) {
        test2 <- test2[!grepl(remove[rm], test2$ssciname,ignore.case = TRUE),] 
      }
      if (nrow(test2) > 1) {test <- test2}
    }
    max <- max(test$pident)
    upper <- as.numeric(test$pident.min.best[1]) # EES set the minimum similarity threshold for including a hit to the minimum similarity of the best matching taxid
    lower <- max-as.numeric(lower_margin)
    test <- test[which(test$pident >= lower),] # select all lines for a query
    test$margin <- "lower"
    test[test$pident >= upper,"margin"] <- "upper"
    
    new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
    i=i+1
  }
  return(new_IDtable)
}

#Function2
# Get full taxonomic path for all hits within the upper limit of each ASV. Identical species are only queried once.

get_classification <- function(IDtable2) {
  require(taxizedb)
  all_staxids <- names(table(IDtable2$staxid[IDtable2$margin=="upper"])) # get all taxids for table
  all_classifications <- list() # prepare list for taxize output
  o=length(all_staxids) # number of taxids
  
  Start_from <- 1 # change if loop needs to be restarted due to time-out
  
  no_taxid_matches <- c()
    
  #Get ncbi classification of each entry
  for (cl in Start_from:o) { # the taxize command "classification" can be run on the all_staxids vector in one line, but often there is
    #a timeout command, therefor this loop workaround.
    print(paste0("step 1 of 3: processing: ", cl , " of ", o , " taxids")) # make a progressline (indicating the index the loops needs to be
    #restarted from if it quits)
     tax_match <- classification(all_staxids[cl], db = "ncbi")   # AGR 
    if (is.na(tax_match) == TRUE) { #This happens when the hit was to a taxon from the BOLD database, which was not in NCBIs taxonomy, when the 
    # BOLD+NCBI database was made, and which were therefore given a dummy taxid. These hits will produce NAs in the classification, or a wrong
    # classification if the dummy taxid is now in use by NCBI. These hits must therefore be manually checked in the output file by comparing the 
    # ssciname column with the final identification (EES 19-08-2022)
      no_taxid_matches <- c(no_taxid_matches,all_staxids[cl])
      tax_na <- classification(all_staxids[1], db = "ncbi") # Use classification of 1st taxid as a template
      tax_na[[1]]["name"]<-"NA" # Replace tax names with NA
      tax_na[[1]]["id"]<-"NA" # Replace taxids at each tax level with NA
      all_classifications[cl] <- tax_na
    } else {
      all_classifications[cl] <- tax_match
    }                                                            
  }
    
  #Construct a taxonomic path from each classification
  output <- data.frame(staxid=character(),kingdom=character(), phylum=character(),class=character(),order=character(),family=character(),genus=character(),species=character(), stringsAsFactors=FALSE)
  totalnames <- length(all_staxids)
  for (curpart in seq(1:totalnames)) {
    print(paste0("step 2 of 3: progress: ", round(((curpart/totalnames) * 100),0) ,"%")) # make a progress line
    currenttaxon <- all_classifications[curpart][[1]]
    if (nchar(currenttaxon[1]) > 0) {
      spec <- all_staxids[curpart]
      output[curpart,"kingdom"] <- currenttaxon[which(currenttaxon$rank == "kingdom"),"name"][1]
      output[curpart,"phylum"] <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"][1]
      output[curpart,"class"] <- currenttaxon[which(currenttaxon$rank == "class"),"name"][1]
      output[curpart,"order"] <- currenttaxon[which(currenttaxon$rank == "order"),"name"][1]
      output[curpart,"family"] <- currenttaxon[which(currenttaxon$rank == "family"),"name"][1]
      output[curpart,"genus"] <- currenttaxon[which(currenttaxon$rank == "genus"),"name"][1]
      output[curpart,"species"] <- currenttaxon[which(currenttaxon$rank == "species"),"name"][1]
      output[curpart,"staxid"] <-  spec # add that row to the filtered IDtable
    }
  }
  taxonomic_info <- merge(IDtable2,output,by = "staxid", all=TRUE)
  taxonomic_info$species[is.na(taxonomic_info$species)] <- taxonomic_info$ssciname[is.na(taxonomic_info$species)] 
  return(list(taxonomic_info,no_taxid_matches))
}

# Function3
# Function for evaluating the taxonomic assignment of each ASV. All hits within the upper margin are used in the evaluation weighted by their evalue, so that suboptimal matches have a lower weight. All hits within the lower margin are put into the output (but not used for evaluating classification)
evaluate_classification <- function(classified) {
  require(tidyr)
  require(dplyr)
  ids <- names(table(classified$qseqid))
  i <- 1
  for (name in ids) {
    print(paste0("last step: progress: ", round(((i/length(ids)) * 100),0) ,"%")) # make a progressline
    test <- classified[which(classified$qseqid == name),]
    test2 <- test %>% filter(margin == "upper")
    test2$score <- 100*(1/test2$evalue)/sum(1/test2$evalue)  # HERE THE SCORE FOR ALL MATCHES PER ASV IS CALCULATED
    test4 <- test2 %>% filter(margin == "upper") %>%
      dplyr::select(margin,qseqid,sseqid,staxid,pident,score,qcovs,kingdom,phylum,class,order,family,genus,species) %>% 
      group_by(qseqid,kingdom, phylum,class,order,family,genus,species) %>% 
      mutate(species_score=sum(score)) %>% 
      group_by(qseqid,kingdom, phylum,class,order,family,genus) %>% 
      mutate(genus_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum,class,order,family) %>% 
      mutate(family_score=sum(score))%>%
      group_by(qseqid,kingdom, phylum,class,order) %>% 
      mutate(order_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum,class) %>% 
      mutate(class_score=sum(score)) %>%
      group_by(qseqid,kingdom, phylum) %>% 
      mutate(phylum_score=sum(score)) %>%
      group_by(qseqid,kingdom) %>% 
      mutate(kingdom_score=sum(score)) %>% ungroup() %>%
      arrange(-kingdom_score,-phylum_score,-class_score,-order_score,-family_score,-genus_score,-species_score)
    test3 <- test4 %>% slice(1)
    test5 <- test4 %>% distinct(qseqid,sseqid,pident,qcovs,kingdom,phylum,class,order,family,genus,species,kingdom_score,phylum_score,class_score,order_score,family_score,genus_score,species_score) 
    string1 <- test %>% dplyr::select(species,pident) %>% 
      distinct(species,pident) %>% arrange(-pident) %>% t()
    string2 <- toString(unlist(string1))
    test3$alternatives <- string2
    if (i == 1){result <- test3} else {
      result <- rbind(result,test3)
    }
    if (i == 1){result2 <- test2} else {
      result2 <- rbind(result2,test2)
    }
    if (i == 1){result3 <- test5} else {
      result3 <- rbind(result3,test5)
    }
    i=i+1
  }
  total_result <- list(taxonon_table = result, all_taxa_table=result2, all_taxa_table_summed=result3)
  return(total_result)
}

# Classify your ASVs by running the wrapper (functionX) "assign_taxonomy" like this. Consider whether it makes sense to remove specific hits or taxa from the evaluation (see explanation below):
my_classified_result <- assign_taxonomy(IDtable, lower_margin = 2, remove = c("uncultured", "environmental")) # EES removed constant upper_margin specification

tax_table <- my_classified_result$classified_table 

# Add maximum similarity for the best matching taxid to tax_table from summary table (EES 25-01-2022)
tax_table$pident.max.best<-"NA"
for (i in unique (tax_table$qseqid)){
   tax_table[tax_table$qseqid==i,]$pident.max.best<-summary[summary$qseqid==i,]$pident.max.best[1]
}

# Add indication of possible misidentification to tax_table from summary table (EES 25-08-2022)
tax_table$possible.misid<-"NA"
for (i in unique (tax_table$qseqid)){
   tax_table[tax_table$qseqid==i,]$possible.misid<-summary[summary$qseqid==i,]$possible.misid[1]
}
                                      
# Add scientific name of top BLAST hit to tax_table from summary table. This is to inspect hits to taxa in 
# the BOLD database, which were not in the NCBI taxonomy when the BOLD+NCBI database was made, and which
# were therefore given a dummy taxid. These hits will produce NAs in the classification, or a wrong
# classification if the dummy taxid is now in use by NCBI (EES 19-08-2022)
tax_table$ssciname<-"NA"
for (i in unique (tax_table$qseqid)){
   tax_table[tax_table$qseqid==i,]$ssciname<-summary[summary$qseqid==i,]$ssciname[1]
}
                                      
# Determine a "final" taxonomic ID, using scores combined with a minimum similarity threshold of 98% for species-level id
score.id = c()
for(i in tax_table$qseqid){
  idx =  tax_table$qseqid==i
  vec=c( isTRUE(all.equal(tax_table[idx,]$species_score, 100)) & as.numeric(tax_table[idx,]$pident.max.best)>=98,
  isTRUE(all.equal(tax_table[idx,]$genus_score, 100)) & (as.numeric(tax_table[idx,]$pident.max.best)<98 | tax_table[idx,]$species_score<100),
  isTRUE(all.equal(tax_table[idx,]$family_score, 100)) & tax_table[idx,]$genus_score<100,
  isTRUE(all.equal(tax_table[idx,]$order_score, 100)) & tax_table[idx,]$family_score<100,
  isTRUE(all.equal(tax_table[idx,]$class_score, 100)) & tax_table[idx,]$order_score<100,
  isTRUE(all.equal(tax_table[idx,]$phylum_score, 100)) & tax_table[idx,]$class_score<100)

  # If the species score is 100 and the sequence similarity above 98%, final ID will be to species level
  # If a species level ID cannot be made, but the genus score is 100, final ID will be to genus level
  # If a genus level ID cannot be made, but the family score is 100, final ID will be to family level
  # If a family level ID cannot be made, but the order score is 100, final ID will be to order level
  # If an order level ID cannot be made, but the class score is 100, final ID will be to class level
  # If a class level ID cannot be made, but the phylum score is 100, final ID will be to phylum level
  # If a phylum level ID cannot be made, but the kingdom score is 100, final ID will be to kingdom level

  condition_num = which(vec)[1]
  columns = c("species","genus","family","order","class","phylum")

  if(is.na(condition_num))
    score.id = c(score.id, "NA")
  else if(condition_num==1)
    score.id = c(score.id, as.vector(tax_table[idx,]$species)[1])
  else if(condition_num==2)
    score.id = c(score.id, as.vector(tax_table[idx,]$genus)[1])
  else if(condition_num==3)
    score.id = c(score.id, as.vector(tax_table[idx,]$family)[1])
  else if(condition_num==4)
    score.id = c(score.id, as.vector(tax_table[idx,]$order)[1])
  else if(condition_num==5)
    score.id = c(score.id, as.vector(tax_table[idx,]$class)[1])
  else if(condition_num==6)
    score.id = c(score.id, as.vector(tax_table[idx,]$phylum)[1])
  else
    score.id = c(score.id, as.vector(tax_table[idx,]$kingdom)[1])
}

tax_table$score.id <- score.id
   
# Optionally, synonyms of scientific names can be downloaded from an appropriate database (WoRMS in the example below. See "taxize" documentation for other database options).
# However, it should be stressed that such a search should be complemented with manual searching across databases, as it is unlikely to be exhaustive.
      
# Get World Register of Marine Species (WoRMS) IDs and search the WoRMS database for synonyms

# Detach taxizedb package before loading taxize, as these packages have similar functions     
#detach("package:taxizedb", unload=TRUE)

# Load taxize package
#library(taxize) 

# Create a column providing the number of words in the final identification of each sequence.     
#tax_table$nwords <- sapply(strsplit(as.character(tax_table$score.id), " "), length)
      
# Create a column for the synonyms found in the chosen database    
#tax_table$synonyms<-"NA"  

# Create a column for the currently valid scientific name
#tax_table$valid_name<-"NA"

# For each eDNA sequence with a species-level identification (nwords equal to 2), add synonyms and valid name     
#for (i in tax_table$qseqid){
#
#   if (tax_table[idx,]$nwords==2){
#       
#       worms_id<-get_wormsid(tax_table[idx,]$species)
#
#       syns<-synonyms(as.character(worms_id[1]),db="worms")
#       
#       if (dim(syns[[as.character(worms_id[1])]])>0) {
#       
#                                                       syns_vec<-syns[[as.character(worms_id[1])]]['scientificname']
#
#                                                       tax_table[idx,]$synonyms<-paste(unlist(syns_vec), collapse = ', ')
#                                                       
#                                                       tax_table[idx,]$valid_name<-syns[[as.character(worms_id[1])]][['valid_name']][1]
#                                                       }
#       }
#}

#Write the taxonomic results to a table
write.table(tax_table, args[3], sep = "\t", quote = F, row.names = F)
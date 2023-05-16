#### 1. Load packages and import data ####

install.packages("ggplot2")
library(ggplot2)
install.packages("plyr")
library("plyr")
install.packages("corrplot")
library("corrplot")
install.packages("gllvm")
library("gllvm")
library(gllvm)
install.packages("lifecycle")
library(lifecycle)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
a
library("phyloseq")



#load data - phyloseq
COSQ_rare2 <- readRDS("/work/769436/COSQ_rare2_correct.rds")
COSQ_rare2

#load data - environmental data
sed_metadata <- read.delim("/work/701023/sed_metadata.txt")
wat_metadata <- read.delim("/work/701023/wat_metadata.txt")



#tjekker at formatet af phyloseq objekterne ser rigtigt ud
COSQ_rare2@otu_table[1,c(1:5)] #taxa = rows, colums=samples
COSQ_rare2@tax_table[1,c(1:5)] #taxa = rows, colums = taxonomy
#View(COSQ_rare2@sam_data) #samples = rows, colums = abiotic variables



#### 2. Filtering and optimizing data before modelling ####

###Create extra variables in sample_data
COSQ_rare2@sam_data$sshc <- paste(COSQ_rare2@sam_data$substrate_type,
                                  COSQ_rare2@sam_data$season,
                                  COSQ_rare2@sam_data$habitat,
                                  COSQ_rare2@sam_data$cluster, sep="_") #ssch = substrate, season, cluster, habitat

COSQ_rare2@sam_data$ssc<-paste(COSQ_rare2@sam_data$substrate_type,
                               COSQ_rare2@sam_data$season, 
                               COSQ_rare2@sam_data$cluster, sep="_") #ssc = substrate, season, cluster

COSQ_rare2@sam_data$stc<-paste(COSQ_rare2@sam_data$substrate_type, 
                               COSQ_rare2@sam_data$cluster, sep="_") #stc = substrate type, cluster

COSQ_rare2@sam_data$snch<-paste(COSQ_rare2@sam_data$season, 
                                COSQ_rare2@sam_data$cluster, 
                                COSQ_rare2@sam_data$habitat,sep="_") #snch = season, cluster, habitat


#Subset by substrate type & remove cluster 2 (cluster 2 was only sampled in 1 season)
COSQ_rare2_no_c2 <- subset_samples(COSQ_rare2,!cluster==2) 
doi<-data.frame(sample_data(COSQ_rare2_no_c2)) #new sample data without cluster 2

#COMPARE 
COSQ_rare2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6183 taxa and 936 samples ]
#sample_data() Sample Data:       [ 936 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6183 taxa by 8 taxonomic ranks ]

COSQ_rare2_no_c2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6183 taxa and 930 samples ]
#sample_data() Sample Data:       [ 930 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6183 taxa by 8 taxonomic ranks ]



### Import and optimize water data (OBS this only contains data from the water samples, not bottom samples)

#Temperature, salinity mv: 
#wat_metadata <- read.delim("C:/Users/Karoline/OneDrive/Speciale/D. Databehandling NYT data/TSV_CSV_TEXT_etc/Tilsendt/wat_metadata.txt")
#View(wat_metadata)


# The snch variable produced above is used to add the salinity and chlorophyll data to the sample data (metadata) created earlier
# Creates a new column "Salinity" in doi (metadata) and imports the salinity column from merged_data to doi
#based on identifcal values in the snch column in the two dataframes (season, habitat, cluster)
#Same procedure for chlorophyll

#add new variables from the water data to the doi sample data
doi$Salinity<-wat_metadata$Salinity[match(doi$snch, wat_metadata$snch)]
doi$Temperature<-wat_metadata$Temperature[match(doi$snch, wat_metadata$snch)]
doi$PO4<-wat_metadata$PO4[match(doi$snch, wat_metadata$snch)]
doi$NO3<-wat_metadata$PO4[match(doi$snch, wat_metadata$snch)]
doi$Chlorophyll<-wat_metadata$Chlorophyll[match(doi$snch, wat_metadata$snch)]
for(i in 1:nrow(doi)){doi[i,"PO4"]<-ifelse(doi[i,"PO4"]<=0,0,doi[i,"PO4"])} #checks whether any values in the PO4 column are less than or equal to zero --> if yes it sets the value to zero. Otherwise it leaves the value unchanged.



# SEDIMENT DATA
#sed_metadata <- read.delim("C:/Users/Karoline/OneDrive/Speciale/D. Databehandling NYT data/TSV_CSV_TEXT_etc/Tilsendt/sed_metadata.txt")
#View(sed_metadata)

#calculate mean value of organic content for each sample (snch = season, cluster, habitat == sample)
OC_f_wat <- ddply(sed_metadata, .(snch), summarise, grp.mean=mean(Organic_content))

#Insert the organic content data from the sediment data into doi (sample data total)
#if the substrate type = sediment --> insert organic_content data from sed_metadata into doi sample data
#if substrate_type is water --> insert the mean-value of organic content for the given sample (the values created in OC_f_wat)
doi$Organic_content<-ifelse(doi$substrate_type=="sediment", sed_metadata$Organic_content[match(doi$root, sed_metadata$Sample_ID)], OC_f_wat$grp.mean[match(doi$snch, OC_f_wat$snch)]) 

#round to five decimals
doi$Organic_content<-round(doi$Organic_content, digits=5)

#mean value of total phospor pr sample
TP_f_wat <- ddply(sed_metadata, .(snch), summarise, grp.mean=mean(TP)) #ddply function --> combine results into a dataframe (for subset of a dataframe)

#Insert TP (total phosphor) into the doi data - same procedure as above for the organic content
doi$TP<-ifelse(doi$substrate_type=="sediment", sed_metadata$TP[match(doi$root, sed_metadata$Sample_ID)], TP_f_wat$grp.mean[match(doi$snch, TP_f_wat$snch)])

doi$TP<-round(doi$TP, digits=4)



#Import water content data to the doi sample data
WC_f_wat <- ddply(sed_metadata, .(snch), summarise, grp.mean=mean(Water_content))

doi$Water_content<-ifelse(doi$substrate_type=="sediment", sed_metadata$Water_content[match(doi$root, sed_metadata$Sample_ID)], WC_f_wat$grp.mean[match(doi$snch, WC_f_wat$snch)])

doi$Water_content<-round(doi$Water_content, digits=4)

#Import organic content
IC_f_wat <- ddply(sed_metadata, .(snch), summarise, grp.mean=mean(Inorganic_content))

doi$Inorganic_content<-ifelse(doi$substrate_type=="sediment", sed_metadata$Inorganic_content[match(doi$root, sed_metadata$Sample_ID)], IC_f_wat$grp.mean[match(doi$snch, IC_f_wat$snch)])

doi$Inorganic_content<-round(doi$Inorganic_content, digits=5)







##Update OTU table
#Identify samples that have no salinity data
with_NA1<-rownames(doi[is.na(doi$Salinity),])
#no samples lack salinity data


#Identify samples that have no chlorophyll data (there are 6)
with_NA2<-rownames(doi[is.na(doi$Chlorophyll),])
#autumn_18_sand and autumn_33_sand lack chlorophyll data


#Remove samples that have missing data for either salinity or chlorophyll - creates a new phyloseq object
newddw<-subset_samples(COSQ_rare2_no_c2, !(root %in% c(with_NA1, with_NA2))) #obs I changed sample_root in Marcelos code to snch as the with_NA2 vector contains season_cluster_habitat names in stead of root names with the new wat_metadata
newddw
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6183 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6183 taxa by 8 taxonomic ranks ]

#Remove OTUs that are no longer present in any samples
tudao0 = filter_taxa(newddw, function(x) sum(x) > 0, TRUE) 
tudao0
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6128 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6128 taxa by 8 taxonomic ranks ]

#EES: Remove empty samples - this seems unnecessary? Check if sample number is the same before and after
#checked: sample number is the same = this step is unnecessary
fdt = prune_samples(sample_sums(tudao0)>0,tudao0)

#CHECK phyloseq object for number of taxa in OTU table
fdt
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6128 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6128 taxa by 8 taxonomic ranks ]

#Agglomerate (merge) taxa from the same phylum
po_dataF<-tax_glom(fdt, "phylum") #reads pr sample in each phylum are summed
po_dataF 
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 42 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 42 taxa by 8 taxonomic ranks ]

#In the original datas�t (pident70) we had 47 phyla
levels(as.factor(pident70_marine$phylum))
#[1] "Annelida"          "Apicomplexa"       "Apusozoa"          "Arthropoda"       
#[5] "Ascomycota"        "Bacillariophyta"   "Basidiomycota"     "Bigyra"           
#[9] "Bryozoa"           "Cercozoa"          "Chlorophyta"       "Choanozoa"        
#[13] "Chordata"          "Chytridiomycota"   "Cnidaria"          "Cryptophyta"      
#[17] "Ctenophora"        "Discosea"          "Echinodermata"     "Endomyxa"         
#[21] "Entoprocta"        "Evosea"            "Gastrotricha"      "Gnathostomulida"  
#[25] "Haptophyta"        "Kinorhyncha"       "Loukozoa"          "Mollusca"         
#[29] "Mucoromycota"      "NA"                "Nematoda"          "Nemertea"         
#[33] "Ochrophyta"        "Onychophora"       "Oomycota"          "Palpitia"         
#[37] "Placozoa"          "Platyhelminthes"   "Porifera"          "Prasinodermophyta"
#[41] "Rhodophyta"        "Rotifera"          "Streptophyta"      "Tardigrada"       
#[45] "Tubulinea"         "Xenacoelomorpha"   "Zoopagomycota"    

#in the Imported COSQ file where we have filtered out non-marine groups and NA's in phylum
#we have 42 taxa
levels(as.factor(COSQ_rare2@tax_table[,2]))
#[1] "Annelida"          "Apicomplexa"       "Apusozoa"          "Arthropoda"       
#[5] "Ascomycota"        "Bacillariophyta"   "Basidiomycota"     "Bigyra"           
#[9] "Bryozoa"           "Cercozoa"          "Chlorophyta"       "Choanozoa"        
#[13] "Chordata"          "Cnidaria"          "Cryptophyta"       "Ctenophora"       
#[17] "Discosea"          "Echinodermata"     "Endomyxa"          "Entoprocta"       
#[21] "Evosea"            "Gastrotricha"      "Gnathostomulida"   "Haptophyta"       
#[25] "Loukozoa"          "Mollusca"          "Mucoromycota"      "Nematoda"         
#[29] "Nemertea"          "Ochrophyta"        "Oomycota"          "Palpitia"         
#[33] "Placozoa"          "Platyhelminthes"   "Porifera"          "Prasinodermophyta"
#[37] "Rhodophyta"        "Rotifera"          "Streptophyta"      "Tardigrada"       
#[41] "Tubulinea"         "Xenacoelomorpha"  

#phyla removed from pident70
#"Chytridiomycota", "Kinorhyncha", "NA", "Onychophora", "Zoopagomycota" 

#in the po_dataF where we have filtered out sediment samples and samples lacking cl and salinity data
#we have 41 taxa
levels(as.factor(po_dataF@tax_table[,2]))
#[1] "Annelida"          "Apicomplexa"       "Apusozoa"          "Arthropoda"       
#[5] "Ascomycota"        "Bacillariophyta"   "Basidiomycota"     "Bigyra"           
#[9] "Bryozoa"           "Cercozoa"          "Chlorophyta"       "Choanozoa"        
#[13] "Chordata"          "Cnidaria"          "Cryptophyta"       "Ctenophora"       
#[17] "Discosea"          "Echinodermata"     "Endomyxa"          "Entoprocta"       
#[21] "Evosea"            "Gastrotricha"      "Gnathostomulida"   "Haptophyta"       
#[25] "Loukozoa"          "Mollusca"          "Mucoromycota"      "Nematoda"         
#[29] "Nemertea"          "Ochrophyta"        "Oomycota"          "Placozoa"         
#[33] "Platyhelminthes"   "Porifera"          "Prasinodermophyta" "Rhodophyta"       
#[37] "Rotifera"          "Streptophyta"      "Tardigrada"        "Tubulinea"        
#[41] "Xenacoelomorpha"  

#phyla in the original phyloseq object COSQ_rare2 that has been removed in po_dataF:
#Palpitia (Chromista)


#Richness
otuo <- data.frame(otu_table(fdt))
otuo1 <- t(data.frame(otuo, check.names=F)) #transpose table, so columns = taxa instead of samples
taxo<-data.frame(tax_table(fdt), stringsAsFactors=FALSE)
colnames(otuo1)<-taxo$phylum #columns in the OTU , OBS changed Phylum to phylum
do<-data.frame(sample_data(fdt))
clades<-levels(factor(taxo$phylum)) #object with all the different phylas in the taxonomy table
clades
otuo2<-t(data.frame(otuo1, check.names=F)) #transpose otuo, columns are now samples again
tabr<-cbind(taxo, otuo2) 
tabr <- tabr[,-3:-6] #OBS, I remove different row numbers compared to Marcelos script, as taxa have different row numbers in ours 
tabr <- tabr[,-4] #remove score.id column
#kept only the kingdom, phylum and species columns fra taxa kolonnen


ch<-do$root #extraxts the column with the root names from ,OBS changes sample_root to root
z<-expand.grid("root"=ch, stringsAsFactors=FALSE) #returns a data frame with the ch vector in it



for(i in 1:length(clades))
{
  gtr<-subset(tabr, phylum==clades[i]) #subsets the phylum column in the tabr dataframe, looks for all the elements (phylas) in the clades object
  rich<-colSums(gtr[,-1:-3] != 0) #sums the numbers in the columns in gtr (sum of reads in each samples) - it calculates the total number of non-zero values in each column, removes column no. 1-3 with taxonomy!OBS I changed this code from Marcelos script
  z<-cbind(z, rich) #puts all elements from the rich object into the z vector with all the sample roots as rows
  t<-1+i #number of phyla + 1 (39+1=40) - what is this for?
  colnames(z)[t]<-clades[i] #inserts the phylum names in clades as column names in the z dataframe
}
rich_asv<-z[,-1] #removes the first column from z, which is "root"
#rich_asv is a df where columns = phyla, rows = samples
#each cell represents the sum of reads in each sample for the given phylum
rich_asv1 <- rich_asv #for later comparison with the other rich_asv object further down in the script

#Saves as table
#write.table(rich_asv1,"rich_asv1_gllvm_script",sep="\t",row.names=F)


#Rel_abund (relative abundance)
datarg = transform_sample_counts(po_dataF, function(x) x/sum(x)) #calculates the relative abundance of reads in each sample - i think it is the reads in each sample pr phylum is divided with the sum of all reads pr sample. In the phyloseq object po_dataF (the input in this code), is where all taxa are aggregated in phylas. The otu table now contains relative abundances in each sample instead of reads pr phylum in each sample.
#View(datarg@otu_table)
datarg
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 42 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 42 taxa by 8 taxonomic ranks ]
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
#View(otu) #rows = taxa
#otu<-t(data.frame(otu,check.names=F)) I remove this code from Marcelos script, as rows and columns are already at the right position in our data
tab<-cbind(tax, otu)
rownames(tab)<-tab$phylum 
tab<-tab[,-1:-8] #remove all taxonomy columns
ttab<-t(data.frame(tab, check.names=F)) #transpose - columns are now = phyla


#Filtering by rich, abund and sites occurring
pabund<-as.data.frame(cbind(log10(colMeans(ttab)), colnames(ttab))) #Percent abundance: ttab is the df with the relative abundance pr sample for each phylum, this code calculates the mean relative abundance in ALL samples pr phylum. Log10 is used to transform numbers to decimal numbers (without e)?- log10 used to bring our values on to a comparable scale
colnames(pabund)[2]<-"clades_p"
prich<-as.data.frame(cbind(log10(colMeans(rich_asv)), colnames(rich_asv))) #richness: The rich_asv df contains the sum of reads pr sample pr phylum.This code calculates the mean value of the reads pr sample pr phyla (mean value of the columns), and binds it with the column names which are the phylum names. Why are the values different than in pabund?
colnames(prich)[2]<-"clades_p"
trich<-t(rich_asv) #transposes rich_asv
prich$sqrt_sites_occur <-sqrt(rowSums(trich != 0)) #sums up the no. of reads pr phylum in all samples and takes the square_root - inserts as a new column in prich. Square root transformation can be usefl for normalizing a skewed distribution. Brings values into a linear relations-ship
abund_rich_summary <- as.data.frame(merge(prich, pabund, by="clades_p")) #creates a table with abundance, richness and sites occuring for each phylum
colnames(abund_rich_summary)<-c("clades_p","log10_mean_rich","sqrt_sites_occur","log10_mean_rel_abund") #new column names
abund_rich_summary$log10_mean_rich<-as.numeric(abund_rich_summary$log10_mean_rich) #mean richness values are transformed to numeric values 
abund_rich_summary$log10_mean_rel_abund<-as.numeric(abund_rich_summary$log10_mean_rel_abund) #relative mean abundance values are transformed to numeric values
abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),] #mean relative abundance values are sorted from lowest to highest values, OBS the values are negative
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),] #same for mean richness
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),] #same for the square root of sites occuring

#saves as table
write.table(abundsorted,"701023/abundsorted_gllvm",sep="\t",row.names=F)


list_top15<-c(abundsorted$clades_p[1:15],richsorted$clades_p[1:15],sitessorted$clades_p[1:15]) #the first 15 rows are filtered out - the 15 phyla with the lowest values of mean relative abundance, richness and occurrence pr site - the 15 phyla with lowest values from each of the listed dataframe = 45 phyla names in total
#View(list_top15)
phy_to_remove<-as.character(unique(list_top15)) #merges all rows that contains the same phylum name = now we have 15 phylas
phy_to_remove #list with the lowest abundant, rich, occuring phylas in our data. 16 in total
#[1] "Entoprocta"        "Endomyxa"          "Palpitia"          "Choanozoa"        
#[5] "Tardigrada"        "Tubulinea"         "Apicomplexa"       "Bigyra"           
#[9] "Cercozoa"          "Bryozoa"           "Xenacoelomorpha"   "Loukozoa"         
#[13] "Mucoromycota"      "Streptophyta"      "Prasinodermophyta" "Gnathostomulida"  
#[17] "Ctenophora"        "Nemertea"   

#add phyla to remove list
#plusP<-c("Ctenophora") #OBS I commented out this Ctenophora as Marcelo recommended us not to remove this phylum unless it causes us troubles

#po_data1F = subset_taxa(po_dataF, !(phylum %in% c(phy_to_remove, plusP)))
po_data1F = subset_taxa(po_dataF, !(phylum %in% c(phy_to_remove))) #I removed the plusP object
po_data1F #there are now 24 phyla in this phyloseq object, as we filtered out the 15 lowest abundant phyla 
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 24 taxa by 8 taxonomic ranks ]


#Filter phyla with abundance lower than 10 in at least 10 % of samples
#Add this step to class modelling
po_data2 = filter_taxa(po_data1F, function(x) sum(x > 10) > (0.1*length(x)), TRUE) #x=reads pr sample (for each phylum. phylums = rows). so the sum function sums up all values for each column = sample. step1: 0.1*length(x) calculates 10 % of the samples (x), length(x) = the number of samples. Step 2: sum(x > 10) takes out all phylas with a sum of reads higher than 10. Together: sum(x > 10) > (0.1*length(x)) - filters out/removes phylas with abundance lower than 10 in at least 10 % of samples
po_data2 #we lose one phyla now - Placozoa has been removed
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 23 taxa by 8 taxonomic ranks ]

#Keep only top 30 phyla (richness)
toprichsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"], decreasing=T ),] #order the data - highest to lowest values of mean richness (decreasing=T)
list_top30<-c(toprichsorted$clades_p[1:30]) #takes out the 30 phyla with the highest richness
list_top30
#[1] "Oomycota"        "Ochrophyta"      "Bacillariophyta" "Discosea"        "Cnidaria"       
#[6] "Arthropoda"      "Rhodophyta"      "Porifera"        "Mollusca"        "Annelida"       
#[11] "Basidiomycota"   "Chlorophyta"     "Rotifera"        "Placozoa"        "Echinodermata"  
#[16] "Cryptophyta"     "Ascomycota"      "Nematoda"        "Chordata"        "Apusozoa"       
#[21] "Haptophyta"      "Evosea"          "Platyhelminthes" "Gastrotricha"    "Mucoromycota"   
#[26] "Bryozoa"         "Bigyra"          "Loukozoa"        "Xenacoelomorpha" "Nemertea"      
po_data2.1 = subset_taxa(po_data2, (phylum %in% list_top30)) #creates a new phyloseq object of the po_data2, that contains only phyla in the top30 list, to make sure, there are no phylas in it, that are not in the top30 richest
po_data2.1 #however, we lose no phyla by this - still 22
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 23 taxa by 8 taxonomic ranks ]

tudao1 = filter_taxa(po_data2, function(x) sum(x) > 0, TRUE) #filters all taxa(phyla) with a sum of samples higher than zero
o_dataF = prune_samples(sample_sums(tudao1)>0,tudao1) #remove all samples with a samples, with a sample sum of zero
o_dataF #no samples has been removed by "prune_samples" as there have been no samples with a sample sum of zero in tudao1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 23 taxa by 8 taxonomic ranks ]


taxoF<-data.frame(tax_table(o_dataF), stringsAsFactors=FALSE) #creates a dataframe of the taxonomy table in the o_dataF
phy_to_keep<-taxoF$phylum #hiver alle taxa (phyla) ud fra o_dataF objektet ovenfor, som er blevet renset for samples og phyla
phy_to_keep #contains the names of all the 23 phyla to keep
#[1] "Bacillariophyta" "Cryptophyta"     "Haptophyta"      "Ochrophyta"      "Oomycota"       
#[6] "Ascomycota"      "Basidiomycota"   "Annelida"        "Arthropoda"      "Chordata"       
#[11] "Cnidaria"        "Echinodermata"   "Mollusca"        "Nematoda"        "Placozoa"       
#[16] "Platyhelminthes" "Porifera"        "Rotifera"        "Chlorophyta"     "Rhodophyta"     
#[21] "Apusozoa"        "Discosea"        "Evosea"      

po_data2 = subset_taxa(fdt, phylum %in% phy_to_keep) #takes the fdt phyloseq object (where empty samples and OTU's are removed, line 140), and filters only the phyla that are contained in the phy_to_keep list
po_data2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 5972 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 5972 taxa by 8 taxonomic ranks ]

tudao1 = filter_taxa(po_data2, function(x) sum(x) > 0, TRUE) #creates a new tudao1 object based on the po_data2 object in the code above, filters all taxa with a sample sum higher than zero
o_data = prune_samples(sample_sums(tudao1)>0,tudao1) #new o_data object, all empty samples are removed
o_data
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 5972 taxa and 918 samples ]
#sample_data() Sample Data:       [ 918 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 5972 taxa by 8 taxonomic ranks ]

#OBS this is the same code as previous in the script - we run this code on the filtered data now. 
otuo<-data.frame(otu_table(o_data)) #new otuo object based on the o_data above
otuo_t <- t(data.frame(otuo, check.names=F)) #I added this code to convert our otu table to the same format as Marcelos. Rows are samples and columns are taxa (phylum)
taxo<-data.frame(tax_table(o_data), stringsAsFactors=FALSE) #new taxo object
colnames(otuo_t)<-taxo$phylum 
ncol(otuo_t) #2739 colums in otuo_t = 2739 taxa
do<-data.frame(sample_data(o_data))


#add water variables to sample data
do$Salinity<-doi$Salinity[match(do$snch, doi$snch)]
do$NO3<-doi$NO3[match(do$snch, doi$snch)]
do$PO4<-doi$PO4[match(do$snch, doi$snch)]
do$Temperature<-doi$Temperature[match(do$snch, doi$snch)]
do$Chlorophyll<-doi$Chlorophyll[match(do$snch, doi$snch)]

#add sediment variables to sample data
do$Organic_content<-doi$Organic_content[match(do$root, doi$root)]
do$TP<-doi$TP[match(do$root, doi$root)]
do$Water_content<-doi$Water_content[match(do$root, doi$root)]
do$Inorganic_content<-doi$Inorganic_content[match(do$root, doi$root)]

do$habitat <- factor(do$habitat, levels = c("sand", "rocks", "eelgrass")) #adds habitat column to sample data
do$season <- factor(do$season, levels = c("spring", "autumn")) #adds season column to sample data
do$cluster<-as.character(do$cluster) #adds cluster column to sample data
do[,12:20]<-scale(do[,12:20]) #column 12-20 are numeric values - scale/normalize


#Mk dummy variables. EES: Convert categorical variables to binary to allow modelling with gllvm
head(do)
do$season<-as.character(do$season)
do$substrate_type<-as.character(do$substrate_type)

for (h in 1:nrow(do))
{
  if(do[h,"season"]=="spring") {
    do[h,"season"]<-as.numeric(0)
  } else {
    do[h,"season"]<-as.numeric(1)}
  
  if(do[h,"substrate_type"]=="sediment") {
    do[h,"substrate_type"]<-as.numeric(0)
  } else {
    do[h,"substrate_type"]<-as.numeric(1)}
  
  if(do[h,"habitat"]=="rocks") {
    do[h,"habitat_rocks"]<-as.numeric(1)
  } else {
    do[h,"habitat_rocks"]<-as.numeric(0)}
  
  if(do[h,"habitat"]=="eelgrass") {
    do[h,"habitat_eel"]<-as.numeric(1)
  } else {
    do[h,"habitat_eel"]<-as.numeric(0)}
  
}

#autumn = 1, spring = 0
#habitat_rocks = 1 if there is 

do$season<-as.numeric(do$season)
do$substrate_type<-as.numeric(do$substrate_type)

#import date to the do sample data
do$Time<-wat_metadata$Time[match(do$snch, wat_metadata$snch)]

#convert dates to numeric values to the number of days after"31-12-18"
reference_date <- as.Date("31-12-18",format = "%d-%m-%y")
do$dates <- as.Date(do$Time,format = "%d-%m-%y")
do$days_since_reference_date <- as.numeric(do$dates-reference_date)

do$days_since_reference_date <- as.numeric(do$days_since_reference_date)
#check that the 1st of jan has no. 1
jan1 <- as.Date("01-01-19",format = "%d-%m-%y")
as.numeric(jan1-reference_date)
#1


clades<-levels(factor(taxo$phylum)) #creates a new clades object with the names of all phylas in the new taxo df (from o_data phyloseq object). There are 20 phyla in clades.
otuo2<-t(data.frame(otuo_t, check.names=F)) #transposes otuo_t, now taxa = rows again
tabr<-cbind(taxo, otuo2) #binds the new taxo and the new otuo, OBS the code is the same as in line 190-202, but the new taxo and otuo2 objects are based on the o_data phyloseq object from line 312 instead of the fdt phyloseq object from line 136

tabr<-tabr[,-1] #remove kingdom
tabr<-tabr[,-2:-7] #remove all taxonomy except phylum
ch<-do$root #creates a vector containing all 441 sample names
z<-expand.grid("root"=ch, stringsAsFactors=FALSE) #creates a dataframe of ch


for(i in 1:length(clades)) #clades contains all the 22 phyla names also present in tabr (in the phylum column), length=22
{
  gtr<-subset(tabr, phylum==clades[i]) #finds all the phylum names from clades (one at a time) in the tabr dataframe
  rich<-colSums(gtr[,-1] != 0) #sums the no. of reads in all sample for each phylum - removes column 1 = phylum name. Removes all samples phyla with a sample sum of 0?
  z<-cbind(z, rich) #binds the z vector containing the phylum names, and the rich df containg the sum of reads pr sample
  t<-1+i #what is this for? t=40 - maybe for inserting the right column names (one at the time) in the code below
  colnames(z)[t]<-clades[i] #inserts the phylum names from the clades object as column names in the z dataframe
}

rich_asv<-z[,-1]

#write.table(rich_asv,"rich_asv_both_substrates",sep="\t",row.names=F)


#### 3. GLLVM Models - find optimal distribution and no. of lv####

#Test which distribution is more suitable (two latent var)

fit_p2 <- gllvm(rich_asv, do, num.lv = 2,
                formula = ~ Temperature + substrate_type + 
                  habitat_eel + habitat_rocks,
                row.eff = ~ (1 | cluster), 
                family = poisson(), 
                method = "VA", 
                control.start = list(starting.val ="zero"))
fit_p2
#log-likelihood:  -43638.9 
#Residual degrees of freedom:  20953 
#AIC:  87599.8 
#AICc:  87602.29 
#BIC:  88880.99 

saveRDS(fit_p2,"701023/fit_p2.rds")


fit_NB2 <- gllvm(rich_asv, do, num.lv = 2,
                 formula = ~ Temperature + substrate_type + 
                   habitat_eel + habitat_rocks,
                 row.eff = ~ (1 | cluster), 
                 family = "negative.binomial", 
                 method = "VA", 
                 control.start = list(starting.val ="zero"))
fit_NB2
#log-likelihood:  -44638.76 
#Residual degrees of freedom:  20930 
#AIC:  89645.51 
#AICc:  89648.77 
#BIC:  91109.73 

saveRDS(fit_NB2,"701023/fit_NB2.rds")

#Looks like poisson is the most suitable distribution --> gives the lowest AIC value

#Find the model with the optimal number of latent variables - from tutorial

fit_p1 <- gllvm(rich_asv, do, num.lv = 1,
                formula = ~ Temperature + substrate_type + 
                  habitat_eel + habitat_rocks,
                row.eff = ~ (1 | cluster), 
                family = poisson(), 
                method = "VA", 
                control.start = list(starting.val ="zero"))
fit_p1
#log-likelihood:  -49130.76 
#Residual degrees of freedom:  20975 
#AIC:  98539.52 
#AICc:  98541.38 
#BIC:  99645.64 

saveRDS(fit_p1,"701023/fit_p1.rds")

fit_p3 <- gllvm(rich_asv, do, num.lv = 3,
                formula = ~ Temperature + substrate_type + 
                  habitat_eel + habitat_rocks,
                row.eff = ~ (1 | cluster), 
                family = poisson(), 
                method = "VA", 
                control.start = list(starting.val ="zero"))
fit_p3
#log-likelihood:  -42581.28 
#Residual degrees of freedom:  20932 
#AIC:  85526.56 
#AICc:  85529.74 
#BIC:  86974.86 

#Warning message:
# In gllvm(rich_asv, do, num.lv = 3, formula = ~Temperature + substrate_type +  :
#           The algorithm did not converge, the maximum number of iterations might have been reached.

saveRDS(fit_p3,"701023/fit_p3.rds")


#CONCLUSION: The optimal number of latent variables is 3 and the most suitable distribution is poission

#### 4. GLLVM models - Including interactions ####

fit_p3_int <- gllvm(rich_asv, do, num.lv = 3,
                    formula = ~ Temperature*substrate_type*
                      habitat_eel*habitat_rocks,
                    row.eff = ~ (1 | cluster), 
                    family = poisson(), 
                    method = "VA", 
                    control.start = list(starting.val ="zero"))
fit_p3_int
#log-likelihood:  -42596.72 
#Residual degrees of freedom:  20679 
#AIC:  86063.45 
#AICc:  86081.79 
#BIC:  89525.04 

#Warning message:
# In gllvm(rich_asv, do, num.lv = 3, formula = ~Temperature * substrate_type *  :
#           Determinant of the variance-covariance matix is zero. Please double check your model for e.g. overfitting or lack of convergence.
#saveRDS(fit_p3_int,"701023/fit_p3_int.rds")

fit_p3_int1 <- gllvm(rich_asv, do, num.lv = 3,
                     formula = ~ Temperature+substrate_type+
                       habitat_eel+habitat_rocks+
                       Temperature*(substrate_type+habitat_eel+habitat_rocks)+
                       substrate_type*(habitat_eel+habitat_rocks),
                     row.eff = ~ (1 | cluster), 
                     family = poisson(), 
                     method = "VA", 
                     control.start = list(starting.val ="zero"))
fit_p3_int1
#log-likelihood:  -42556.44 
#Residual degrees of freedom:  20817 
#AIC:  85706.87 
#AICc:  85715.38 
#BIC:  88070.31 

#saveRDS(fit_p3_int1,"701023/fit_p3_int1.rds")

#Higher AIC indicate that including interactions does not improve the model.


#### 5. GLLVM - Testing the explanatory variables ####

###Removing one ex. var. at the time

#rem. Temperature
fit_p3_noTemp <- gllvm(rich_asv, do, num.lv = 3,
                       formula = ~ substrate_type + habitat_eel + habitat_rocks,
                       row.eff = ~ (1 | cluster), 
                       family = poisson(), 
                       method = "VA", 
                       control.start = list(starting.val ="zero"))

fit_p3_noTemp
#log-likelihood:  -38572.12 
#Residual degrees of freedom:  20955 
#AIC:  77462.24 
#AICc:  77464.67 
#BIC:  78727.52 

#saveRDS(fit_p3_noTemp,"701023/fit_p3_noTemp.rds")

#rem. substrate type
fit_p3_nosub <- gllvm(rich_asv, do, num.lv = 3,
                      formula = ~ Temperature + habitat_eel + habitat_rocks,
                      row.eff = ~ (1 | cluster), 
                      family = poisson(), 
                      method = "VA", 
                      control.start = list(starting.val ="zero"))

fit_p3_nosub
#log-likelihood:  -42577.05 
#Residual degrees of freedom:  20955 
#AIC:  85472.1 
#AICc:  85474.53 
#BIC:  86737.38

#Warning message:
# In gllvm(rich_asv, do, num.lv = 3, formula = ~Temperature + habitat_eel +  :
#           The algorithm did not converge, the maximum number of iterations might have been reached.

#saveRDS(fit_p3_nosub,"701023/fit_p3_nosub.rds")

#rem. both habitat variables
fit_p3_nohab <- gllvm(rich_asv, do, num.lv = 3,
                      formula = ~ Temperature + substrate_type,
                      row.eff = ~ (1 | cluster), 
                      family = poisson(), 
                      method = "VA", 
                      control.start = list(starting.val ="zero"))

fit_p3_nohab
#log-likelihood:  -42490.01 
#Residual degrees of freedom:  20978 
#AIC:  85252.03 
#AICc:  85253.81 
#BIC:  86334.27 

#saveRDS(fit_p3_nohab,"701023/fit_p3_nohab.rds")

#rem. ellgrass
fit_p3_noeel <- gllvm(rich_asv, do, num.lv = 3,
                      formula = ~ Temperature + substrate_type + habitat_rocks,
                      row.eff = ~ (1 | cluster), 
                      family = poisson(), 
                      method = "VA", 
                      control.start = list(starting.val ="zero"))

fit_p3_noeel
#log-likelihood:  -42428.98 
#Residual degrees of freedom:  20955 
#AIC:  85175.96 
#AICc:  85178.39 
#BIC:  86441.23 

#saveRDS(fit_p3_noeel,"701023/fit_p3_noeel.rds")

#rem. rock
fit_p3_norock <- gllvm(rich_asv, do, num.lv = 3,
                       formula = ~ Temperature + substrate_type + habitat_eel,
                       row.eff = ~ (1 | cluster), 
                       family = poisson(), 
                       method = "VA", 
                       control.start = list(starting.val ="zero"))
fit_p3_norock
#log-likelihood:  -42533.14 
#Residual degrees of freedom:  20955 
#AIC:  85384.27 
#AICc:  85386.7 
#BIC:  86649.54 

#Warning message:
# In gllvm(rich_asv, do, num.lv = 3, formula = ~Temperature + substrate_type +  :
#           The algorithm did not converge, the maximum number of iterations might have been reached.

#saveRDS(fit_p3_norock,"701023/fit_p3_norock.rds")


###Inluding only 1 ex. var.

#Temp
fit_p3_Temp <- gllvm(rich_asv, do, num.lv = 3,
                     formula = ~ Temperature,
                     row.eff = ~ (1 | cluster), 
                     family = poisson(), 
                     method = "VA", 
                     control.start = list(starting.val ="zero"))

fit_p3_Temp
#log-likelihood:  -42663.22 
#Residual degrees of freedom:  21001 
#AIC:  85552.44 
#AICc:  85553.66 
#BIC:  86451.66 
#saveRDS(fit_p3_Temp,"701023/fit_p3_Temp.rds")


#substrate_type
fit_p3_sub <- gllvm(rich_asv, do, num.lv = 3,
                    formula = ~ substrate_type,
                    row.eff = ~ (1 | cluster), 
                    family = poisson(), 
                    method = "VA", 
                    control.start = list(starting.val ="zero"))

fit_p3_sub
#log-likelihood:  -38933.19 
#Residual degrees of freedom:  21001 
#AIC:  78092.38 
#AICc:  78093.61 
#BIC:  78991.6 
#saveRDS(fit_p3_sub,"701023/fit_p3_sub.rds")


#habitat
fit_p3_hab <- gllvm(rich_asv, do, num.lv = 3,
                    formula = ~ habitat_eel + habitat_rocks,
                    row.eff = ~ (1 | cluster), 
                    family = poisson(), 
                    method = "VA", 
                    control.start = list(starting.val ="zero"))

fit_p3_hab
#log-likelihood:  -39868.12 
#Residual degrees of freedom:  20978 
#AIC:  80008.24 
#AICc:  80010.02 
#BIC:  81090.49 
#saveRDS(fit_p3_hab,"701023/fit_p3_hab.rds")

#habitat_eel
fit_p3_eel <- gllvm(rich_asv, do, num.lv = 3,
                    formula = ~ habitat_eel,
                    row.eff = ~ (1 | cluster), 
                    family = poisson(), 
                    method = "VA", 
                    control.start = list(starting.val ="zero"))

fit_p3_eel
#log-likelihood:  -40083.24 
#Residual degrees of freedom:  21001 
#AIC:  80392.49 
#AICc:  80393.72 
#BIC:  81291.71 
#saveRDS(fit_p3_eel,"701023/fit_p3_eel.rds")


#habitat_rocks
fit_p3_rocks <- gllvm(rich_asv, do, num.lv = 3,
                      formula = ~ habitat_rocks,
                      row.eff = ~ (1 | cluster), 
                      family = poisson(), 
                      method = "VA", 
                      control.start = list(starting.val ="zero"))

fit_p3_rocks
#log-likelihood:  -39955.89 
#Residual degrees of freedom:  21001 
#AIC:  80137.77 
#AICc:  80139 
#BIC:  81036.99 
#saveRDS(fit_p3_rocks ,"701023/fit_p3_rocks.rds")


### Including 2 ex. var

#Temp and sub
fit_p3_TempSub <- gllvm(rich_asv, do, num.lv = 3,
                        formula = ~ Temperature + substrate_type,
                        row.eff = ~ (1 | cluster), 
                        family = poisson(), 
                        method = "VA", 
                        control.start = list(starting.val ="zero"))
fit_p3_TempSub
#log-likelihood:  -42490 
#Residual degrees of freedom:  20978 
#AIC:  85251.99 
#AICc:  85253.77 
#BIC:  86334.24 
#saveRDS(fit_p3_TempSub,"701023/fit_p3_TempSub.rds")


#Temp and habitat_eel
fit_p3_TempEel <- gllvm(rich_asv, do, num.lv = 3,
                        formula = ~ Temperature + habitat_eel,
                        row.eff = ~ (1 | cluster), 
                        family = poisson(), 
                        method = "VA", 
                        control.start = list(starting.val ="zero"))
fit_p3_TempEel
#log-likelihood:  -42818.82 
#Residual degrees of freedom:  20978 
#AIC:  85909.65 
#AICc:  85911.42 
#BIC:  86991.89 
#saveRDS(fit_p3_TempEel,"701023/fit_p3_TempEel.rds")

#Temp and habitat_rocks
fit_p3_TempRocks <- gllvm(rich_asv, do, num.lv = 3,
                          formula = ~ Temperature + habitat_rocks,
                          row.eff = ~ (1 | cluster), 
                          family = poisson(), 
                          method = "VA", 
                          control.start = list(starting.val ="zero"))

fit_p3_TempRocks
#log-likelihood:  -42605.16 
#Residual degrees of freedom:  20978 
#AIC:  85482.32 
#AICc:  85484.1 
#BIC:  86564.57
#saveRDS(fit_p3_TempRocks,"701023/fit_p3_TempRocks.rds")


#substrate type and habitat_eel
fit_p3_SubEel <- gllvm(rich_asv, do, num.lv = 3,
                       formula = ~ substrate_type + habitat_eel,
                       row.eff = ~ (1 | cluster), 
                       family = poisson(), 
                       method = "VA", 
                       control.start = list(starting.val ="zero"))
fit_p3_SubEel
#log-likelihood:  -38821.1 
#Residual degrees of freedom:  20978 
#AIC:  77914.21 
#AICc:  77915.98 
#BIC:  78996.45 
#saveRDS(fit_p3_SubEel,"701023/fit_p3_SubEel.rds")


#substrate type and habitat_rocks

fit_p3_SubRocks <- gllvm(rich_asv, do, num.lv = 3,
                         formula = ~ substrate_type + habitat_rocks,
                         row.eff = ~ (1 | cluster), 
                         family = poisson(), 
                         method = "VA", 
                         control.start = list(starting.val ="zero"))
fit_p3_SubRocks
#log-likelihood:  -38673.2 
#Residual degrees of freedom:  20978 
#AIC:  77618.4 
#AICc:  77620.17 
#BIC:  78700.64 
#saveRDS(fit_p3_SubRocks,"701023/fit_p3_SubRocks.rds")



#habitat_eel and habitat_rocks
fit_p3_hab <- gllvm(rich_asv, do, num.lv = 3,
                    formula = ~ habitat_eel + habitat_rocks,
                    row.eff = ~ (1 | cluster), 
                    family = poisson(), 
                    method = "VA", 
                    control.start = list(starting.val ="zero"))

fit_p3_hab
#log-likelihood:  -39868.12 
#Residual degrees of freedom:  20978 
#AIC:  80008.24 
#AICc:  80010.02 
#BIC:  81090.49 
#saveRDS(fit_p3_hab,"701023/fit_p3_hab.rds")


### No explanatory variables
fit_p3_noex <- gllvm(rich_asv, num.lv = 3,
                     formula = ~1,
                     family = poisson(), 
                     method = "VA", 
                     control.start = list(starting.val ="zero"))
fit_p3_noex
#log-likelihood:  -40461.41 
#Residual degrees of freedom:  21025 
#AIC:  81100.83 
#AICc:  81101.59 
#BIC:  81809.06 
#saveRDS(fit_p3_noex,"701023/fit_p3_noex.rds")



### Season and date included instead of temp

#season
fit_p3_season <- gllvm(rich_asv, do, num.lv = 3,
                       formula = ~ season + substrate_type + 
                         habitat_eel + habitat_rocks,
                       row.eff = ~ (1 | cluster), 
                       family = poisson(), 
                       method = "VA", 
                       control.start = list(starting.val ="zero"))

#Warning message:
#  In gllvm(rich_asv, do, num.lv = 3, formula = ~season + substrate_type +  :
#            The algorithm did not converge, the maximum number of iterations might have been reached.

fit_p3_season
#log-likelihood:  -38471.53 
#Residual degrees of freedom:  20932 
#AIC:  77307.05 
#AICc:  77310.24 
#BIC:  78755.35 

#saveRDS(fit_p3_season,"701023/fit_p3_season.rds")

#only season
fit_p3_onlyseason <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season,
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))
fit_p3_onlyseason
#log-likelihood:  -40052.62 
#Residual degrees of freedom:  21001 
#AIC:  80331.24 
#AICc:  80332.47 
#BIC:  81230.46 

#saveRDS(fit_p3_onlyseason,"701023/fit_p3_onlyseason.rds")


#including date
fit_p3_date <- gllvm(rich_asv, do, num.lv = 3,
                     formula = ~ days_since_reference_date + substrate_type + 
                       habitat_eel + habitat_rocks,
                     row.eff = ~ (1 | cluster), 
                     family = poisson(), 
                     method = "VA", 
                     control.start = list(starting.val ="zero"))
fit_p3_date
#log-likelihood:  -38630.43 
#Residual degrees of freedom:  20932 
#AIC:  77624.87 
#AICc:  77628.05 
#BIC:  79073.17 

#saveRDS(fit_p3_date,"701023/fit_p3_date.rds")


#only date included

fit_p3_onlydate <- gllvm(rich_asv, do, num.lv = 3,
                         formula = ~ days_since_reference_date,
                         row.eff = ~ (1 | cluster), 
                         family = poisson(), 
                         method = "VA", 
                         control.start = list(starting.val ="zero"))
fit_p3_onlydate
#log-likelihood:  -39980.3 
#Residual degrees of freedom:  21001 
#AIC:  80186.59 
#AICc:  80187.82 
#BIC:  81085.81 

#saveRDS(fit_p3_onlydate,"701023/fit_p3_date.rds")


#### 6. season gllvm ####

#season
fit_p3_season <- gllvm(rich_asv, do, num.lv = 3,
                       formula = ~ season + substrate_type + 
                         habitat_eel + habitat_rocks,
                       row.eff = ~ (1 | cluster), 
                       family = poisson(), 
                       method = "VA", 
                       control.start = list(starting.val ="zero"))

#Warning message:
#  In gllvm(rich_asv, do, num.lv = 3, formula = ~season + substrate_type +  :
#            The algorithm did not converge, the maximum number of iterations might have been reached.

fit_p3_season
#log-likelihood:  -38471.53 
#Residual degrees of freedom:  20932 
#AIC:  77307.05 
#AICc:  77310.24 
#BIC:  78755.35 


#interactions
fit_p3_seasonint <- gllvm(rich_asv, do, num.lv = 3,
                          formula = ~ season*substrate_type* 
                            habitat_eel*habitat_rocks,
                          row.eff = ~ (1 | cluster), 
                          family = poisson(), 
                          method = "VA", 
                          control.start = list(starting.val ="zero"))

fit_p3_seasonint
#log-likelihood:  -38418.52 
#Residual degrees of freedom:  20679 
#AIC:  77707.05 
#AICc:  77725.39 
#BIC:  81168.64 
#saveRDS(fit_p3_seasonint,"701023/fit_p3_seasonint.rds")

#AIC increases in the interaction model - we do not include


#only season interactions
fit_p3_seasonint1 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*(substrate_type + 
                                                 habitat_eel + habitat_rocks),
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint1
#log-likelihood:  -38386.64 
#Residual degrees of freedom:  20863 
#AIC:  77275.29 
#AICc:  77281.35 
#BIC:  79272.67 
#saveRDS(fit_p3_seasonint1,"701023/fit_p3_seasonint1.rds")



#season habitat int
fit_p3_seasonint2 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*(habitat_eel + habitat_rocks),
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint2
#log-likelihood:  -39716.73 
#Residual degrees of freedom:  20909 
#AIC:  79843.46 
#AICc:  79847.5 
#BIC:  81474.78 
#saveRDS(fit_p3_seasonint2,"701023/fit_p3_seasonint2.rds")

#season sub rock int
fit_p3_seasonint3 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*(substrate_type + habitat_rocks),
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint3
#log-likelihood:  -38403.17 
#Residual degrees of freedom:  20909 
#AIC:  77216.34 
#AICc:  77220.38 
#BIC:  78847.67
#saveRDS(fit_p3_seasonint3,"701023/fit_p3_seasonint3.rds")


#season sub eel int
fit_p3_seasonint4 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*(substrate_type + habitat_eel),
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint4
#log-likelihood:  -38577.68 
#Residual degrees of freedom:  20909 
#AIC:  77565.35 
#AICc:  77569.39 
#BIC:  79196.68 
#saveRDS(fit_p3_seasonint4,"701023/fit_p3_seasonint4.rds")


#season sub int
fit_p3_seasonint5 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*substrate_type,
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint5
#log-likelihood:  -38679.45 
#Residual degrees of freedom:  20955 
#AIC:  77676.89 
#AICc:  77679.32 
#BIC:  78942.16 
#saveRDS(fit_p3_seasonint5,"701023/fit_p3_seasonint5.rds")

#season eel int
fit_p3_seasonint6 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*habitat_eel,
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint6
#log-likelihood:  -39978.43 
#Residual degrees of freedom:  20955 
#AIC:  80274.86 
#AICc:  80277.29 
#BIC:  81540.13 
#saveRDS(fit_p3_seasonint6,"701023/fit_p3_seasonint6.rds")


#season rock int
fit_p3_seasonint7 <- gllvm(rich_asv, do, num.lv = 3,
                           formula = ~ season*habitat_rocks,
                           row.eff = ~ (1 | cluster), 
                           family = poisson(), 
                           method = "VA", 
                           control.start = list(starting.val ="zero"))

fit_p3_seasonint7
#log-likelihood:  -39793.52 
#Residual degrees of freedom:  20955 
#AIC:  79905.03 
#AICc:  79907.46 
#BIC:  81170.31 
#saveRDS(fit_p3_seasonint7,"701023/fit_p3_seasonint7.rds")


#season + hab
fit_p3_seasonhab <- gllvm(rich_asv, do, num.lv = 3,
                          formula = ~ season + 
                            habitat_eel + habitat_rocks,
                          row.eff = ~ (1 | cluster), 
                          family = poisson(), 
                          method = "VA", 
                          control.start = list(starting.val ="zero"))
fit_p3_seasonhab
#log-likelihood:  -39758.74 
#Residual degrees of freedom:  20955 
#AIC:  79835.48 
#AICc:  79837.9 
#BIC:  81100.75
#saveRDS(fit_p3_seasonhab,"701023/fit_p3_seasonhab.rds")

#season + hab_eel
fit_p3_seasonhab_eel <- gllvm(rich_asv, do, num.lv = 3,
                              formula = ~ season + 
                                habitat_eel,
                              row.eff = ~ (1 | cluster), 
                              family = poisson(), 
                              method = "VA", 
                              control.start = list(starting.val ="zero"))
fit_p3_seasonhab_eel
#log-likelihood:  -40832.2 
#Residual degrees of freedom:  20978 
#AIC:  81936.41 
#AICc:  81938.19 
#BIC:  83018.66 
#saveRDS(fit_p3_seasonhab_eel,"701023/fit_p3_seasonhab_eel.rds")

#season + hab_rocks
fit_p3_seasonhab_rocks <- gllvm(rich_asv, do, num.lv = 3,
                                formula = ~ season + 
                                  habitat_rocks,
                                row.eff = ~ (1 | cluster), 
                                family = poisson(), 
                                method = "VA", 
                                control.start = list(starting.val ="zero"))
fit_p3_seasonhab_rocks
#log-likelihood:  -39811.64 
#Residual degrees of freedom:  20978 
#AIC:  79895.28 
#AICc:  79897.06 
#BIC:  80977.52 
#saveRDS(fit_p3_seasonhab_rocks,"701023/fit_p3_seasonhab_rocks.rds")


#season + substratetype
fit_p3_seasonsub <- gllvm(rich_asv, do, num.lv = 3,
                          formula = ~ season + substrate_type,
                          row.eff = ~ (1 | cluster), 
                          family = poisson(), 
                          method = "VA", 
                          control.start = list(starting.val ="zero"))
fit_p3_seasonsub
#log-likelihood:  -38848.6 
#Residual degrees of freedom:  20978 
#AIC:  77969.19 
#AICc:  77970.97 
#BIC:  79051.44 
#saveRDS(fit_p3_seasonsub,"701023/fit_p3_seasonsub.rds")


#### 7. date gllvm ####
fit_p3_date <- gllvm(rich_asv, do, num.lv = 3,
                     formula = ~ days_since_reference_date + substrate_type + 
                       habitat_eel + habitat_rocks,
                     row.eff = ~ (1 | cluster), 
                     family = poisson(), 
                     method = "VA", 
                     control.start = list(starting.val ="zero"))
fit_p3_date
#log-likelihood:  -38630.43 
#Residual degrees of freedom:  20932 
#AIC:  77624.87 
#AICc:  77628.05 
#BIC:  79073.17 

#saveRDS(fit_p3_date,"701023/fit_p3_date.rds")

#date interaction model
fit_p3_dateint <- gllvm(rich_asv, do, num.lv = 3,
                        formula = ~ days_since_reference_date*substrate_type* 
                          habitat_eel*habitat_rocks,
                        row.eff = ~ (1 | cluster), 
                        family = poisson(), 
                        method = "VA", 
                        control.start = list(starting.val ="zero"))

fit_p3_dateint
#log-likelihood:  -38440.52 
#Residual degrees of freedom:  20679 
#AIC:  77751.04 
#AICc:  77769.39 
#BIC:  81212.64 
#saveRDS(fit_p3_dateint,"701023/fit_p3_dateint.rds")

#The AIC values increases for the interaction model
#therefore we do not include interactions


#date + hab
fit_p3_datehab <- gllvm(rich_asv, do, num.lv = 3,
                        formula = ~ days_since_reference_date + 
                          habitat_eel + habitat_rocks,
                        row.eff = ~ (1 | cluster), 
                        family = poisson(), 
                        method = "VA", 
                        control.start = list(starting.val ="zero"))
fit_p3_datehab
#log-likelihood:  -40898.76 
#Residual degrees of freedom:  20955 
#AIC:  82115.51 
#AICc:  82117.94 
#BIC:  83380.79
#saveRDS(fit_p3_datehab,"701023/fit_p3_datehab.rds")

#date + hab_eel
fit_p3_datehab_eel <- gllvm(rich_asv, do, num.lv = 3,
                            formula = ~ days_since_reference_date + 
                              habitat_eel,
                            row.eff = ~ (1 | cluster), 
                            family = poisson(), 
                            method = "VA", 
                            control.start = list(starting.val ="zero"))
fit_p3_datehab_eel
#log-likelihood:  -40000.76 
#Residual degrees of freedom:  20978 
#AIC:  80273.52 
#AICc:  80275.3 
#BIC:  81355.77
#saveRDS(fit_p3_datehab_eel,"701023/fit_p3_datehab_eel.rds")

#date + hab_rocks
fit_p3_datehab_rocks <- gllvm(rich_asv, do, num.lv = 3,
                              formula = ~ days_since_reference_date + 
                                habitat_rocks,
                              row.eff = ~ (1 | cluster), 
                              family = poisson(), 
                              method = "VA", 
                              control.start = list(starting.val ="zero"))
fit_p3_datehab_rocks
#log-likelihood:  -39808.87 
#Residual degrees of freedom:  20978 
#AIC:  79889.75 
#AICc:  79891.52 
#BIC:  80971.99 
#saveRDS(fit_p3_datehab_rocks,"701023/fit_p3_datehab_rocks.rds")


#date + substratetype
fit_p3_datesub <- gllvm(rich_asv, do, num.lv = 3,
                        formula = ~ days_since_reference_date + substrate_type,
                        row.eff = ~ (1 | cluster), 
                        family = poisson(), 
                        method = "VA", 
                        control.start = list(starting.val ="zero"))
fit_p3_datesub
#log-likelihood:  -38860.79 
#Residual degrees of freedom:  20978 
#AIC:  77993.58 
#AICc:  77995.35 
#BIC:  79075.82 
#saveRDS(fit_p3_datesub,"701023/fit_p3_datesub.rds")


#### Importing RDS files locally ####
fit_p3_noTemp <- readRDS("C:/Users/Karoline/OneDrive/Speciale/G. Databehandling NYT KORREKT data/RDS filer/gllvm modeller both substrates/fit_p3_noTemp.rds")
fit_p3_noTemp
#log-likelihood:  -38572.12 
#Residual degrees of freedom:  20955 
#AIC:  77462.24 
#AICc:  77464.67 
#BIC:  78228.97 


rescov0 <- getResidualCov(fit_p3_noTemp)
rescov0$trace
rescov0$var.q

plot(fit_p3_noTemp)
dev.off() #clears plot window


#New ugly (but correct) plots

rbPal <- c("#D2B48C", "#5F9EA0")
# "#8B4513", "#D2691E", "#006400", "#32CD32","#B266FF", "#FF66FF")
hbPal <- c("#6666FF","#FF8000","#00CC00")
sbPal2 <- c("#66FFFF","#66CC00","#CCCC00","#FF9933","#FF0000")
cbPal <- c("#B2FF66","#80FF00","#00CC00","#009900","#336600")

#Inspect colour breaks in continuous variables
jajas<-doi$Salinity[match(do$snch, doi$snch)]
binx<-as.numeric(cut(jajas, breaks = 5))


dp<-do
dp$legend<-NA
dp$legend<-ifelse(dp$season==0,as.numeric(0),as.numeric(1))
dp$pchh<-ifelse(dp$habitat=="sand",as.numeric(0),ifelse(dp$habitat=="rocks",as.numeric(1),as.numeric(2)))

pdf("ordi_phylum_nb_wat_rich_null3r_LVs12_ug.pdf")
par(mfrow = c(2,2), mar=c(2,2,2,2))
Col <- rbPal[as.numeric(cut(dp[,"legend"], breaks = 2))]
ordiplot(fit_p3_noTemp, symbols = T, s.colors = Col, main = "season", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("spring","autumn"), col=rbPal, pch=1, bty = "n")
Cols <- sbPal2[as.numeric(cut(dp[,"Salinity"], breaks = 5))]
ordiplot(fit_p3_noTemp, symbols = T, s.colors = Cols, main = "salinity", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("7-11","12-17","18-22","23-27","28-33"), col=sbPal2, pch=1, bty = "n")

Colh <- hbPal[as.numeric(cut(dp[,"pchh"], breaks = 3))]
ordiplot(fit_p3_noTemp, symbols = T, s.colors = Colh, main = "habitat", cex.spp=0.4, which.lvs = c(1,2))
legend("topleft", legend = c("sand", "rocks", "eelgrass"), col=hbPal, pch=1, bty = "n")
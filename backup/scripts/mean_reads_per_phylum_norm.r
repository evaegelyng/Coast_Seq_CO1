#MEAN READS PR PHYLYM NORMALISERET DATA

#### 1. Load data and packages ####
library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")

COSQ_rare2 <- readRDS("C:/Users/au620760/OneDrive - Aarhus universitet/COSQ_results_2023/Karoline/R scripts/COSQ_rare2.rds") #s�dan indl�ser du RDS filen
COSQ_rare2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3207 taxa and 932 samples ]
#sample_data() Sample Data:       [ 932 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 3207 taxa by 8 taxonomic ranks ]


#### 2. Merging OTU and taxonomy table  ####

tax_otu<-merge(COSQ_rare2@tax_table,COSQ_rare2@otu_table,by="row.names")
#OBS Nu kan du se objektet som en dataframe


#### 3. Aggregate phyla ####
n <- ncol(tax_otu)
n #941
phyla_samples <- tax_otu[,c(2,3,10:n)] #hiver kingdom og phylum kolonnen ud samt alle samples kolonnerne
agg_phyla_samples<- aggregate(.~phylum+kingdom,data=phyla_samples,sum) #merger alle r�kker med ens phylum og kingdomnavne -summerer alle counts i sample kolonnerne

#bytter om p� phylum og kingdom kolonnen
agg_phyla_samples <- agg_phyla_samples[,c(2,1,3:934)]


#Jeg laver en ny kolonne "mean counts" som indeholder middelv�rdien af antal counts pr sample
agg_phyla_samples$mean_counts <- rowMeans(agg_phyla_samples[,c(3:934)])
agg_phyla_samples$mean_counts1<- round(agg_phyla_samples$mean_counts) #runder op til helt tal
agg_phyla_samples$mean_counts2<- round(agg_phyla_samples$mean_counts,2) #med 2 decimaler


#Jeg laver ogs� en kolonne med summen af reads pr sample
agg_phyla_samples$sum_counts <- rowSums(agg_phyla_samples[,c(3:934)])

#nyt objekt kun med kingdom, phylum, mean_counts og sum_counts
#1: Kingdom, 2: phylum, 936=mean_counts1, 937=sum_counts, 938 = mean_counts (med to decimaler)
counts_pr_phyla <- agg_phyla_samples[,c(1,2,936,937,938)]

#nye navne til kolonner
names(counts_pr_phyla)[3] <- "mean_reads_pr_sample"
names(counts_pr_phyla)[4] <- "mean_reads_pr_sample_2decimals"
names(counts_pr_phyla)[5] <- "sum_reads_pr_sample"

#ordner i  r�kkef�lge
counts_pr_phyla<- counts_pr_phyla[order(-(counts_pr_phyla$mean_reads_pr_sample)),]
View(counts_pr_phyla)

#Gemmer som tabel
write.table(counts_pr_phyla,"total and mean reads pr sample pr phylum (data normalized twice by median depth)",sep="\t",row.names=F)



# Barplots

library(ggplot2)
bp_sumreads <- ggplot(counts_pr_phyla,aes(reorder(x=factor(phylum),sum_reads_pr_sample), y=sum_reads_pr_sample, fill= factor(kingdom), colour = factor(kingdom)))+
  geom_bar(width=0.8,stat="identity",position=position_dodge())+
  geom_point()+
  theme(axis.title.x = element_text(margin = margin(r = 20)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.3)) +
  labs(y='Sum of reads',x='Phylum',legend='Kingdom')+
  ggtitle("Sum of reads in all samples pr. phylum") 
bp_sumreads


#top 10 phyla plot
top10 <- counts_pr_phyla[1:10,]

bp_top10_sumreads <- ggplot(top10,aes(reorder(x=factor(phylum),sum_counts), y=sum_counts, fill= factor(kingdom), colour = factor(kingdom)))+
  geom_bar(width=0.8,stat="identity",position=position_dodge())+
  geom_point()+
  theme(axis.title.x = element_text(margin = margin(r = 20)))+
  theme(axis.text.x = element_text(angle = 90,vjust=0.3)) +
  labs(y='Sum of reads',x='Phylum',legend='Kingdom')+
  ggtitle("Sum of reads in all samples - Top 10 phyla") 
bp_top10_sumreads

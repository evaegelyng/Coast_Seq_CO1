# Run from MJOLNIR/results/
# conda activate hmsc

library("phyloseq")
library("ggplot2")

#import data with check.names=F
otu_tab<-read.table("COSQ_otu_phyloseq.txt", sep="\t", header=T, row.names=1,check.names=F)
otu_mat<-as.matrix(otu_tab)

#import classification data with check.names=F
mjolnir_output<-read.table("COSQ_final_dataset.tsv", sep="\t", header=T, row.names=1,check.names=F,colClasses="character")
#Extract the classification columns and other non-sample columns. 
n<-ncol(mjolnir_output)
tax_tab <- mjolnir_output[ c(1:26,n) ]
tax_mat<-as.matrix(tax_tab)

meta_tab <- read.table("metadata/COSQ_metadata_complete.txt", sep="\t", header = T)

# Create phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_mat)
SAM = sample_data(meta_tab)

p_SILVA = phyloseq(OTU, TAX_b, SAM)
p_SILVA

# Remove field samples to include only controls
p_SILVA_2<-subset_samples(p_SILVA, !source=="Field_sample")
p_SILVA_2

# Merge PCR reps and field reps
p_SILVA_3<-merge_samples(p_SILVA_2, "sample_root")
p_SILVA_3

## Rebuild sample data, as the merge_samples function only handles merging of the OTU table
d<-data.frame(sample_data(p_SILVA_3)[,c("sample_root","source")])

d$root<-rownames(d)
d$source<-meta_tab$source[match(d$root, meta_tab$sample_root)]

sample_data(p_SILVA_3)<-d[,c("root","source")]

# Merge ASVs by phylum
p_SILVA_4<-tax_glom(p_SILVA_3, taxrank="phylum")
p_SILVA_4

# Make barplot of raw read abundance
barplot<-plot_bar(p_SILVA_4, "source", fill="phylum")+geom_bar(stat="identity")+
  theme_classic() + 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size=14), axis.text = element_text(size = 12),axis.title.x = element_blank(), axis.title.y = element_text(size=14),plot.margin = unit(c(0,0.4,0,0.4), 'lines'))+
  ggtitle("")+
    labs(y = "Read count")

ggsave(barplot,filename="barplot_contaminants.png",width=10,height=10)

# Remove low-abundance classes
TopPhyla <- names(sort(taxa_sums(p_SILVA_4), TRUE)[1:6])
p_SILVA_5 <- prune_taxa(TopPhyla, p_SILVA_4)

# Make barplot of raw read abundance
barplot<-plot_bar(p_SILVA_5, "source", fill="phylum")+geom_bar(stat="identity")+
  scale_fill_manual(breaks = c("Arthropoda","Ascomycota","Chordata"), values = c("#1F78B4","purple","hotpink")) + 
  theme_classic() + 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size=14), axis.text = element_text(size = 12),axis.title.x = element_blank(), axis.title.y = element_text(size=14),plot.margin = unit(c(0,0.4,0,0.4), 'lines'))+
  ggtitle("")+
    labs(y = "Read count")

ggsave(barplot,filename="barplot_contaminants_top6.png",width=10,height=10)
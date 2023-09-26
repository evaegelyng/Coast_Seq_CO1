# First, create environment with the needed packages
# conda create -n hmsc bioconductor-phyloseq r-ggplot2 r-hmsc r-vegan r-reshape2 r-plyr r-scales r-stringr r-RColorBrewer r-corrplot
# Activate environment: conda activate hmsc

library("phyloseq")
library("ggplot2")
library("Hmsc")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
library("corrplot")
#sessionInfo()

# Load rarefied 70%-similarity dataset (MOTU dataset)
COSQ_rare2<-readRDS("results/COSQ_rare2_correct.rds")

# Load metadata
metadatas<-sample_data(COSQ_rare2)
# Create extra metadata variables
metadatas$sshc<-paste(metadatas$substrate_type, metadatas$season, metadatas$habitat, metadatas$cluster, sep="_")
metadatas$ssc<-paste(metadatas$substrate_type, metadatas$season, metadatas$cluster, sep="_")
metadatas$stc<-paste(metadatas$substrate_type, metadatas$cluster, sep="_")
metadatas$ssh<-paste(metadatas$substrate_type, metadatas$season, metadatas$habitat, sep="_")

sampledatas = sample_data(data.frame(metadatas, row.names=metadatas$root, stringsAsFactors=FALSE))

# Make new phyloseq object with extended metadata
OTU_COI = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
p_COI = phyloseq(OTU_COI, TAX_S)
COSQ_new = merge_phyloseq(p_COI, sampledatas)
COSQ_new

# Remove cluster 2, which was only sampled in spring
COSQ_no_c2<-subset_samples(COSQ_new, !cluster==2)
COSQ_no_c2 # Six samples were removed

#Subset to water substrate
COSQ_w<-subset_samples(COSQ_no_c2, substrate_type=="water")

#Load env data
pc_bs<-read.table("results/metadata/merged_metadata_230427.txt", sep="\t", header=T, row.names=1)
fwat<-read.table("results/metadata/wat_metadata.txt", sep="\t", header=T)

#Load distance matrix
distsea<-read.table("results/metadata/dist_by_sea.txt", sep="\t", header=T)
#mds<-acast(distsea, SiteA ~ SiteB)

#####Third test: structure - merged field replicates
DT1<-merge_samples(COSQ_w, "sshc", fun = mean)
DT1

#Fix sample_data
d<-data.frame(sample_data(DT1)[,c("cluster","season","habitat")])
d$habitat<- sapply(strsplit(as.character(rownames(d)), "_"), function(x) unlist(strsplit(x[3], "_"))[1])
d$season<- sapply(strsplit(as.character(rownames(d)), "_"), function(x) unlist(strsplit(x[2], "_"))[1])
d$sshc<-rownames(d)
d$sch<-paste(d$season,d$cluster,d$habitat,sep="_")
d$coast_shape<-ifelse(d$cluster==2|d$cluster==8|d$cluster==9|d$cluster==10|d$cluster==12|d$cluster==16|d$cluster==29|d$cluster==24,"fjord","open")
sample_data(DT1)<-d[,c("sshc","sch","cluster","season","habitat","coast_shape")]

# Normalize by rarefying
#DT1.1<-rarefy_even_depth(DT1,sample.size = min(sample_sums(DT1)), 
#        trimOTUs = TRUE, verbose = TRUE, replace=FALSE, rngseed= 13072021)
# This removed 4524 OTUs! Even when using median depth as rarefaction threshold,
# more than half the OTUs are lost. Perhaps it is sufficient that the data has 
# been normalized at the field replicate level?

#merge at phylum level. NB! Class was used for 16S! But class has NA values here
DT1.2<-tax_glom(DT1, taxrank="phylum")

#check match between sample data and sp. data
#sdr<-rownames(pc_bs)
#sdo<-sample_names(DT1.2)
#sdr==sdo # Do not match, but these lines were commented out by MPA

#Random subsampling taxa and prepare species data
#set.seed(123)
#randomord20 <- sample(taxa_names(DT1.2), 30, replace=FALSE)
#DT1.3 <- prune_taxa(randomord20, DT1.2)

dataYabund<-data.frame(otu_table(DT1.2))
taxa<-data.frame(tax_table(DT1.2))
colnames(dataYabund)<-taxa$phylum

# Subset phyloseq object to the taxa remaining after normalizing
tax_to_keep<-colnames(dataYabund)
DT1.4 = subset_taxa(DT1, (phylum %in% tax_to_keep))
DT1.4 # No change compared to DT1, as we have not normalized

# Prepare for calculating richness of each phylum in each sample
otuo<-data.frame(otu_table(DT1.4))
taxo<-data.frame(tax_table(DT1.4), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$phylum
do<-data.frame(sample_data(DT1.4))
clades<-levels(factor(taxo$phylum))
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)
# Check tax label columns and remove unnecessary ones
tabr[1:2,1:9]
tabr<-within(tabr,rm("kingdom","class","order","family","genus","species","score.id.pident.70"))
ch<-do$sshc
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

for(i in 1:length(clades))
{
gtr<-subset(tabr, phylum==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

rich_asv<-z[,-1]
dataYrich<-rich_asv

# Calculate relative read abundance
datarg = transform_sample_counts(DT1.2, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$phylum
# Check tax label columns and remove unnecessary ones
tab[1:2,1:9]
tab<-within(tab,rm("kingdom","phylum","class","order","family","genus","species","score.id.pident.70"))
ttab<-t(data.frame(tab, check.names=F))
class(ttab) <- "numeric"

######

#Make dist data rectangle
ds<-distsea[with(distsea, order(SiteA, SiteB)), ]
ds2 <- with(ds, Mads_calc) # Mads_calc is not defined anywhere?
ss <- with(ds, unique(c(as.character(SiteA), as.character(SiteB))))
attributes(ds2) <- with(ds, list(Size = length(ss),
                                  Labels = ss,
                                  Diag = FALSE,
                                  Upper = FALSE,
                                  method = "user"))
class(ds2) <- "dist"

#Format date as days from first sampling day
dsds<-sort(as.Date(pc_bs$Time, format="%d-%m-%y"))

pc_bs$Time_d<-as.integer(abs(as.Date(pc_bs$Time, format="%d-%m-%y")-dsds[1]))

#Prepare env data
pc_bs$season <- factor(pc_bs$season, levels = c("spring","autumn"))
pc_bs$habitat <- factor(pc_bs$habitat, levels = c("sand","rocks","eelgrass"))
pc_bs$sc<-paste(pc_bs$season, pc_bs$cluster, sep="_")
pc_bs$sch<-rownames(pc_bs)
pc_bs$sshc<-paste("water",pc_bs$season,pc_bs$habitat, pc_bs$cluster, sep="_")

#Add temperature missing data
#
#Missing data - temperature data imputed from average of adjacent sites
tempdata<-subset(pc_bs, season=="spring"&(cluster==11|cluster==14))

pc_bs["spring_12_rocks","Temperature"]<-mean(tempdata[1:2,"Temperature"])
pc_bs["spring_12_sand","Temperature"]<-mean(tempdata[3:4,"Temperature"])
pc_bs["spring_13_rocks","Temperature"]<-mean(tempdata[1:2,"Temperature"])
pc_bs["spring_13_sand","Temperature"]<-mean(tempdata[3:4,"Temperature"])

pc_bs2<-pc_bs[,c("Salinity","log_Si","log_PO4","log_DN","Temperature","cube_d14N_15N","Water_content","habitat","season","sshc")]

# Identify samples with complete metadata
good_samples<-rownames(pc_bs2)[rowSums(is.na(pc_bs2)) == 0]
pc_bs3<-pc_bs2
rownames(pc_bs3)<-pc_bs3$sshc
good_samples_sp<-rownames(pc_bs3)[rowSums(is.na(pc_bs3)) == 0]

# Subset data to samples with complete metadata
dataYrich <- dataYrich[good_samples_sp, ]
dataYabund <- dataYabund[good_samples_sp, ]

pc_bs2<-na.omit(pc_bs2)

pc_bs2$sch<-rownames(pc_bs2)
Time_d<-data.frame(pc_bs[,c("season","Time_d")])
Time_d$season<-as.numeric(Time_d$season)

##Attention, forcing change in time in order to make unique season_date combinations per cluster
Time_d[4,2]<-131
Time_d[53,2]<-145
Time_d[152,2]<-5

Time_d <- Time_d[good_samples, ]

Time_d2<-as.matrix(unique(Time_d))
btr<-pc_bs
btr<-btr[good_samples, ]
rownames(Time_d2)<-unique(btr$sc)

######    Preparing alternative envdata for analysis with substrate type
dare<-data.frame(sample_data(DT1.2))
dare <- dare[good_samples_sp, ]
dare$sch<-paste(dare$season,dare$cluster, dare$habitat, sep="_")
dare$sc<-paste(dare$season,dare$cluster, sep="_")
dare$habitat <- factor(dare$habitat, levels = c("sand","rocks","eelgrass"))
dare$coast_shape<-factor(dare$coast_shape, levels = c("fjord","open"))
dare$cluster<-as.factor(dare$cluster)
dare$season<-as.factor(dare$season)
dare$sch<-as.factor(dare$sch)
dare$sshc<-as.factor(dare$sshc)
dare$sc<-as.factor(dare$sc)
dare$Salinity<-pc_bs2$Salinity[match(dare$sch, pc_bs2$sch)]
dare$log_Si<-pc_bs2$log_Si[match(dare$sch, pc_bs2$sch)]
dare$log_PO4<-pc_bs2$log_PO4[match(dare$sch, pc_bs2$sch)]
dare$log_DN<-pc_bs2$log_DN[match(dare$sch, pc_bs2$sch)]
dare$Temperature<-pc_bs2$Temperature[match(dare$sch, pc_bs2$sch)]
dare$cube_d14N_15N<-pc_bs2$cube_d14N_15N[match(dare$sch, pc_bs2$sch)]
dare$Water_content<-pc_bs2$Water_content[match(dare$sch, pc_bs2$sch)]

dare<-dare[,c("sch","cluster","season","sc","sshc","habitat","coast_shape",
"Salinity","log_Si","log_PO4","log_DN","Temperature","cube_d14N_15N","Water_content")]

#Filtering by rich, abund and sites occurring
pabund<-as.data.frame(cbind(log10(colMeans(dataYabund)), colnames(dataYabund)))
colnames(pabund)[2]<-"clades_p"
prich<-as.data.frame(cbind(log10(colMeans(dataYrich)), colnames(dataYrich)))
colnames(prich)[2]<-"clades_p"
trich<-t(dataYrich)
prich$sqrt_sites_occur <-sqrt(rowSums(trich != 0))
abund_rich_summary <- as.data.frame(merge(prich, pabund, by="clades_p"))
colnames(abund_rich_summary)<-c("clades_p","log10_mean_rich","sqrt_sites_occur","log10_mean_rel_abund")
abund_rich_summary$log10_mean_rich<-as.numeric(abund_rich_summary$log10_mean_rich)
abund_rich_summary$log10_mean_rel_abund<-as.numeric(abund_rich_summary$log10_mean_rel_abund)
abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

#Get 50% most prevalent phyla
rich_p50<-abund_rich_summary[abund_rich_summary$sqrt_sites_occur > quantile(abund_rich_summary$sqrt_sites_occur, 0.5), ]
tax_to_keep2<-as.character(rich_p50$clades_p)

#Update tables
dataYrich <- dataYrich[, tax_to_keep2]
dataYabund <- dataYabund[, tax_to_keep2]

#Check prevalence and abundance
P = colMeans(dataYrich>0)
A = colSums(dataYrich)

pdf("results/hist_COI_rich_p50prevphylum_water_plus.pdf")
par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "log10 Richness")
dev.off()

P = colMeans(dataYabund>0)
A = colSums(dataYabund)/sum(dataYabund)

pdf("results/hist_abund_COI_p50prevphylum_water_plus.pdf")
par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "log10 Abundance")
dev.off()

studyDesign = data.frame(
    season = as.factor(dare$season),
    sch = as.factor(dare$sch),
    habitat = as.factor(dare$habitat),
    coast_shape = as.factor(dare$coast_shape),
    cluster = as.factor(dare$cluster),
    Time_d = as.factor(dare$sc),
    space = as.factor(dare$cluster),
    sshc = as.factor(dare$sshc)
)

#Create random effects structure
rl1 = HmscRandomLevel(sData = Time_d2)
rl2 = HmscRandomLevel(distMat = ds2)
rlr = HmscRandomLevel(unit = unique(dare$cluster))

rl1 = setPriors(rl1,nfMax=5)
rl2 = setPriors(rl2,nfMax=5)

rnd_ef<-list("Time_d"=rl1,"space"=rl2)

#Create formula
XFormula = ~ poly(Salinity, degree = 2, raw = TRUE) + log_Si + log_PO4 + log_DN + Temperature + 
cube_d14N_15N + Water_content + coast_shape + habitat

#Set models; abund or rich, evaluate different distributions
m.full.r.lp = Hmsc(Y=dataYrich, XData=dare, XFormula=XFormula,
studyDesign=studyDesign, ranLevels=rnd_ef, distr="lognormal poisson")

#Run
nChains = 2
nParallel = nChains
thin = 10
samples = 1000
transient = 0.2*thin*samples
verbose = samples/10

models = sampleMcmc(m.full.r.lp, thin = thin,
samples = samples, transient = transient,
nChains = nChains, nParallel = nParallel, verbose = verbose)

#m.spatial = sampleMcmc(m.spatial, thin = 10, samples = 1000, transient = 1000, nChains = 2, verbose = 0, updater=list(GammaEta=FALSE))
#m.spatial = sampleMcmc(m.spatial, thin = 15, samples = 1500, transient = 1500, nChains = 4, nParallel=4, verbose = 0, updater=list(GammaEta=FALSE))
#mpost = convertToCodaObject(m.spatial)
#summary(mpost$Beta)

saveRDS(models, file = "results/COI_rich_p50prevphylum_water_plus_c2_thin10_s1000_tr02.rds")

#Diagnostic

modelsII<-models

WAIC2 = computeWAIC(hM=modelsII, byColumn=FALSE)

preds = computePredictedValues(modelsII)
MF= evaluateModelFit(hM = modelsII, predY = preds)
mpost = convertToCodaObject(modelsII, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
ess.beta = effectiveSize(mpost$Beta)
gd<-gelman.diag(mpost$Beta, multivariate=T)

WAIC2

#Explanatory power
#For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and
# predicted values, times the sign of the correlation (SR2).
mean(MF$SR2)
sd(MF$SR2)

#Predictive power
mean(ess.beta)
sd(ess.beta)
mean(gd$psrf)
sd(gd$psrf)

pdf("results/diag_COI_rich_p50prevphylum_water_plus_c2_thin10_s1000_tr02.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("results/diag_beta_COI_rich_p50prevphylum_water_plus_c2_thin10_s1000_tr02.pdf")
plot(mpost$Beta)
dev.off()

head(modelsII$X)
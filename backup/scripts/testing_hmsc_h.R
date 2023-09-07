setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc")


#Diagnostic

modelsh2<-readRDS("16S_rich_p50prevclass_water_plus_c2_thin10_s1000_tr02.rds")

models<-modelsh2
modelsII<-models

WAIC2 = computeWAIC(hM=modelsII, byColumn=FALSE)

preds = computePredictedValues(modelsII)
MF= evaluateModelFit(hM = modelsII, predY = preds)
mpost = convertToCodaObject(modelsII, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
ess.beta = effectiveSize(mpost$Beta)
gd<-gelman.diag(mpost$Beta, multivariate=T)

WAIC2

#Explanatory power
mean(MF$SR2)
sd(MF$SR2)

#Predictive power
mean(ess.beta)
sd(ess.beta)
mean(gd$psrf)
sd(gd$psrf)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc/diag")

pdf("diag_rich_p75prevorder_water_plus_c2_thin10_s1000_tr02.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("diag_beta_rich_p75prevorder_water_plus_c2_thin10_s1000_tr02.pdf")
plot(mpost$Beta)
dev.off()

head(modelsII$X)

#Select variables according to substrate type
#water
groupnames = c("salinity","log_Si","log_PO4","log_DN","Temperature","coast_shape","habitat")
group = c(1,1,2,3,4,5,6,7,7)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)

#water_plus
groupnames = c("salinity","log_Si","log_PO4","log_DN","Temperature","cube_d14N_15N","Water_content","coast_shape","habitat")
group = c(1,1,2,3,4,5,6,7,8,9,9)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)

#seediment
groupnames = c("salinity","cube_d14N_15N","log_TP","log_C","Water_content","coast_shape","habitat")
group = c(1,1,2,3,4,5,6,7,7)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)

#sediment plus
groupnames = c("salinity","cube_d14N_15N","Organic_content","Water_content","log_Si","log_DN","habitat")
group = c(1,1,2,3,4,5,6,7,7)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)


#CLASS/ORDER plots
#run phyloseq taxa filtering (Class/Order)

#
taxh = subset(tax, (Order %in% tax_to_keep))
#
taxh = subset(tax, (Class %in% tax_to_keep))

#CAUTION
#If bacteria, fix colnames
colnames(taxh)[1]<-"Division"
#continue..

#
taxh2<-taxh[,c("Division","Phylum","Class","Order")]
#
taxh2<-taxh[,c("Division","Phylum","Class")]

VPtable<-data.frame(t(VPi$vals))

#Import R2
pr2<-cbind(models$spNames,MF$SR2)

#
colnames(pr2)<-c("Order","SR2")
#
colnames(pr2)<-c("Class","SR2")

pr2<-data.frame(pr2)
pr2$SR2<-as.numeric(pr2$SR2)


#CAUTION
#Check colnames
colnames(VPtable)<-c("habitat","salinity","substrate","time","space")

#Alternatively
colnames(VPtable)<-c("habitat","salinity","time","space","substrate")

#
colnames(VPtable)<-c("salinity","cube_d14N_15N","log_TP","log_C",         
"Water_content","coast_shape","habitat","time","space")

#
colnames(VPtable)<-c("salinity","log_Si","log_PO4","log_DN",
"temperature","coast_shape","habitat","time","space")


#
colnames(VPtable)<-c("salinity","log_Si","log_PO4","log_DN",
"temperature","cube_d14N_15N","Water_content","coast_shape","habitat","time","space")

#
colnames(VPtable)<-c("salinity","cube_d14N_15N","Organic_content","Water_content",
"log_Si","log_DN","habitat","time","space")



#
VPtable$Class<-taxh2$Class[match(rownames(VPtable), taxh2$Order)]
VPtable$Phylum<-taxh2$Phylum[match(rownames(VPtable), taxh2$Order)]
VPtable$Division<-taxh2$Division[match(rownames(VPtable), taxh2$Order)]
VPtable$Order<-rownames(VPtable)
VPtable$SR2<-pr2$SR2[match(rownames(VPtable), pr2$Order)]
#
VPtable$Phylum<-taxh2$Phylum[match(rownames(VPtable), taxh2$Class)]
VPtable$Division<-taxh2$Division[match(rownames(VPtable), taxh2$Class)]
VPtable$Class<-rownames(VPtable)
VPtable$SR2<-pr2$SR2[match(rownames(VPtable), pr2$Class)]

#continue
VPtable$residual<-1-VPtable$SR2

#
nVPtable<-VPtable[,1:5]*VPtable[,9]
nVPtable$residual<-VPtable$residual
nVPtable[7:9]<-VPtable[,6:8]
colnames(nVPtable)[7:9]<-c("Phylum","Division","Class")

#Alternatively
nVPtable<-VPtable[,1:9]*VPtable[,14]
nVPtable$residual<-VPtable$residual
nVPtable[11:14]<-VPtable[,10:13]
colnames(nVPtable)[11:14]<-c("Class","Phylum","Division","Order")

#Alternatively
nVPtable<-VPtable[,1:11]*VPtable[,16]
nVPtable$residual<-VPtable$residual
nVPtable[13:16]<-VPtable[,12:15]
colnames(nVPtable)[13:16]<-c("Class","Phylum","Division","Order")

#Alternatively
nVPtable<-VPtable[,1:9]*VPtable[,13]
nVPtable$residual<-VPtable$residual
nVPtable[11:13]<-VPtable[,10:12]
colnames(nVPtable)[11:13]<-c("Phylum","Division","Class")

#Alternatively
nVPtable<-VPtable[,1:9]*VPtable[,13]
nVPtable$residual<-VPtable$residual
nVPtable[11:13]<-VPtable[,10:12]
colnames(nVPtable)[11:13]<-c("Division","Phylum","Class")

#Alternatively
nVPtable<-VPtable[,1:9]*VPtable[,13]
nVPtable$residual<-VPtable$residual
nVPtable[11:13]<-VPtable[,10:12]
colnames(nVPtable)[11:13]<-c("Division","Phylum","Class")

#Alternatively
nVPtable<-VPtable[,1:11]*VPtable[,15]
nVPtable$residual<-VPtable$residual
nVPtable[13:15]<-VPtable[,12:14]
colnames(nVPtable)[13:15]<-c("Phylum","Division","Class")

#Check colnames
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c("habitat","salinity","substrate","time","space","residual"))

#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class","Order"), measure=c("salinity","log_Si","log_PO4","log_DN","temperature","coast_shape",
"habitat","time","space","residual"))

#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class","Order"), measure=c("salinity","log_Si","log_PO4","log_DN","temperature","cube_d14N_15N",
"Water_content","coast_shape","habitat","time","space","residual"))

#or 16S
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c("salinity","log_Si","log_PO4","log_DN","temperature","coast_shape",
"habitat","time","space","residual"))

#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c(
"salinity","habitat","coast_shape","cube_d14N_15N","log_TP","log_C",         
"Water_content","time","space","residual"))

#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class","Order"), measure=c(
"salinity","habitat","coast_shape","cube_d14N_15N","log_TP","log_C",         
"Water_content","time","space","residual"))

#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c("salinity","cube_d14N_15N","Organic_content","Water_content",
"log_Si","log_DN","habitat","time","space","residual"))

#or
#or
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c("salinity","log_Si","log_PO4","log_DN","temperature","cube_d14N_15N",
"Water_content","coast_shape","habitat","time","space","residual"))



#Plot variance partition
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc/tax_VP")

####Select palette colors:
#brewer.pal(10, "RdBu")
#brewer.pal(12, "Paired")
#

VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("residual","habitat","salinity","substrate","time","space"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#999999", "#CC9966", "#006699", "#CC66CC", "#009900", "#666666"),labels = c(
paste("residual",round(mean(nVPtable$residual),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("substrate",round(mean(nVPtable$substrate),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2))))
ggsave("barplot_vp_big_division_p60class_abund_sshc_logp_polysal_thin25_s1000_tr02_I_substraterandom.pdf")


#16S sediment plus
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","cube_d14N_15N","Organic_content","Water_content",
"log_Si","log_DN","habitat","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3",
"#C7EAE5","#80CDC1","#35978F","#01665E","#003C30"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), 
paste("cube_d14N_15N",round(mean(nVPtable$cube_d14N_15N),digits=2)),
paste("Organic_content",round(mean(nVPtable$Organic_content),digits=2)),
paste("Water_content",round(mean(nVPtable$Water_content),digits=2)),
paste("log_Si",round(mean(nVPtable$log_Si),digits=2)),
paste("log_DN",round(mean(nVPtable$log_DN),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)),  paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("16S_rich_p50prevclass_sediment_plus_c2_thin10_s1000_tr02.pdf")

#Or
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class, VPtable2$Order),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$Order<-factor(VPtable2.0$Order, levels=unique(VPtable2.0$Order))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","log_Si","log_PO4","log_DN","temperature","coast_shape",
"habitat","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Order, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Class+Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0",
"#92C5DE","#4393C3","#2166AC","#053061"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("log_Si",round(mean(nVPtable$log_Si),digits=2)), paste("log_PO4",round(mean(nVPtable$log_PO4),digits=2)),
paste("log_DN",round(mean(nVPtable$log_DN),digits=2)),
paste("temperature",round(mean(nVPtable$temperature),digits=2)),
paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("rich_p75prevorder_water_c2_thin10_s1000_tr02.pdf")


#Or
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class, VPtable2$Order),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$Order<-factor(VPtable2.0$Order, levels=unique(VPtable2.0$Order))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","log_Si","log_PO4","log_DN","temperature","cube_d14N_15N",
"Water_content","coast_shape","habitat","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Order, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Class+Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
"#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("log_Si",round(mean(nVPtable$log_Si),digits=2)), paste("log_PO4",round(mean(nVPtable$log_PO4),digits=2)),
paste("log_DN",round(mean(nVPtable$log_DN),digits=2)),
paste("temperature",round(mean(nVPtable$temperature),digits=2)),
paste("cube_d14N_15N",round(mean(nVPtable$cube_d14N_15N),digits=2)),
paste("Water_content",round(mean(nVPtable$Water_content),digits=2)),
paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("rich_p75prevorder_water_plus_c2_thin10_s1000_tr02.pdf")


#Or
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","log_Si","log_PO4","log_DN","temperature","cube_d14N_15N",
"Water_content","coast_shape","habitat","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
"#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("log_Si",round(mean(nVPtable$log_Si),digits=2)), paste("log_PO4",round(mean(nVPtable$log_PO4),digits=2)),
paste("log_DN",round(mean(nVPtable$log_DN),digits=2)),
paste("temperature",round(mean(nVPtable$temperature),digits=2)),
paste("cube_d14N_15N",round(mean(nVPtable$cube_d14N_15N),digits=2)),
paste("Water_content",round(mean(nVPtable$Water_content),digits=2)),
paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("16S_rich_p50prevclass_water_plus_c2_thin10_s1000_tr02.pdf")


#or 16S water
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","log_Si","log_PO4","log_DN","temperature","coast_shape",
"habitat","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#D1E5F0",
"#92C5DE","#4393C3","#2166AC","#053061"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("log_Si",round(mean(nVPtable$log_Si),digits=2)), paste("log_PO4",round(mean(nVPtable$log_PO4),digits=2)),
paste("log_DN",round(mean(nVPtable$log_DN),digits=2)),
paste("temperature",round(mean(nVPtable$temperature),digits=2)),
paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("16S_rich_p50prevorder_water_c2_thin10_s1000_tr02.pdf")


#Or
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class),]
VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","habitat","coast_shape","cube_d14N_15N","log_TP","log_C",         
"Water_content","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3",
"#C7EAE5","#80CDC1","#35978F","#01665E","#003C30"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("cube_d14N_15N",round(mean(nVPtable$cube_d14N_15N),digits=2)),
paste("log_TP",round(mean(nVPtable$log_TP),digits=2)),
paste("log_C",round(mean(nVPtable$log_C),digits=2)),
paste("Water_content",round(mean(nVPtable$Water_content),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("16S_rich_p50prevorder_sediment_c2_thin10_s1000_tr02.pdf")


#Or
VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class, VPtable2$Order),]

VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))
VPtable2.0$Order<-factor(VPtable2.0$Order, levels=unique(VPtable2.0$Order))


VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("salinity","habitat","coast_shape","cube_d14N_15N","log_TP","log_C",         
"Water_content","time","space","residual"))

ggplot(data = VPtable2.0, aes(x = Order, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Class+Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3",
"#C7EAE5","#80CDC1","#35978F","#01665E","#003C30"),labels = c(
paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("coast_shape",round(mean(nVPtable$coast_shape),digits=2)),
paste("cube_d14N_15N",round(mean(nVPtable$cube_d14N_15N),digits=2)),
paste("log_TP",round(mean(nVPtable$log_TP),digits=2)),
paste("log_C",round(mean(nVPtable$log_C),digits=2)),
paste("Water_content",round(mean(nVPtable$Water_content),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2)),
paste("residual",round(mean(nVPtable$residual),digits=2))))
ggsave("rich_p75prevorder_sediment_c2_thin10_s1000_tr02.pdf")




m_pa<-readRDS("abund_pa_p25_order_I.rds")
m_ln_cop<-readRDS("abund_log_cop_p25_order_I.rds")

postBeta_pa = getPostEstimate(m_pa, parName="Beta")
postBeta_ln_cop = getPostEstimate(m_ln_cop, parName="Beta")

head(postBeta_pa$mean)
head(postBeta_ln_cop$mean)

#DiagII
preds<-readRDS("preds_my_test_models_simple_I.rds")
emf<-readRDS("emf_my_test_models_simple_II.rds")


mean(emf$SR2)

postBeta = getPostEstimate(models, parName="Beta")

head(models$X)
head(postBeta$mean)


#plot(beta.RRR[2,], -postBeta$mean[2,])
#abline(0, 1)
preds = computePredictedValues(m_ln_cop)
MF= evaluateModelFit(hM = m_ln_cop, predY = preds)


partition = createPartition(modelsII, nfolds = 3, column = "cluster")

nChains = 3
nParallel = nChains

cvpreds = computePredictedValues(modelsII, partition = partition, nParallel = nParallel, verbose=modelsII$verbose)

preds = computePredictedValues(modelsII)
MF= evaluateModelFit(hM = modelsII, predY = preds)

preds = computePredictedValues(modelsII,
partition = partition)
MFCV= evaluateModelFit(hM = modelsII, predY = preds)


WAIC2_pa = computeWAIC(hM=m_pa, byColumn=FALSE)
WAIC2_ln_cop = computeWAIC(hM=m_ln_cop, byColumn=FALSE)

mpost = convertToCodaObject(modelsII, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

mpost_pa = convertToCodaObject(m_pa, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
mpost_ln_cop = convertToCodaObject(m_ln_cop, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

ess.beta_pa = effectiveSize(mpost_pa$Beta)
ess.beta_ln_cop = effectiveSize(mpost_ln_cop$Beta)

mpost = convertToCodaObject(models)

ess.beta = effectiveSize(mpost$Beta)

ess.alpha = effectiveSize(mpost$Alpha)

ess.rho = effectiveSize(mpost$Rho)

setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc/diag")

pdf("diag_p60class_rich_polysal_thin25_s1000_tr02.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("diag_beta_p60class_rich_polysal_thin25_s1000_tr02.pdf")
plot(mpost$Beta)
dev.off()

pdf("diag_p40class_polysal_thin500_s2000_tr02_prior.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("diag_beta_p40class_polysal_thin500_s2000_tr02_prior.pdf")
plot(mpost$Beta)
dev.off()


pdf("diag_test_p25order_pa.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost_pa$Beta), main="ess(beta)")
hist(gelman.diag(mpost_pa$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("test_diag_beta_p25order_pa.pdf")
plot(mpost_pa$Beta)
dev.off()


pdf("diag_test_p25order_ln_cop.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost_ln_cop$Beta), main="ess(beta)")
hist(gelman.diag(mpost_ln_cop$Beta, multivariate=T)$psrf, main="psrf(beta)")
dev.off()

pdf("test_diag_beta_p25order_ln_cop.pdf")
plot(mpost_ln_cop$Beta)
dev.off()



#OLD
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc")

VP = computeVariancePartitioning(modelsII)


plotVariancePartitioning =
    function (hM, VP, cols=NULL, main = "Variance Partitioning", ...)
{
   ng = dim(VP$vals)[1]
   if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
   }
   leg = VP$groupnames
   for (r in seq_len(hM$nr)) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
   }
   means = round(100 * rowMeans(VP$vals), 1)
   for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
   }
   dal<-data.frame(VP$vals)
   dal2  <- dal[ , order(names(dal))]
   barplot(as.matrix(dal2), main = main, xlab= "", ylab = "Variance proportion", las = 1,
           legend.text = leg, args.legend= list(x = "topleft", cex=0.9), col = cols,...)
#   mtext("Species", 1,line = 1)
}

#Check colnames
head(modelsII$X)

groupnames = c("habitat","salinity","substrate")
group = c(2,2,1,1,3)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)

groupnames = c("habitat","salinity","substrate")
group = c(2,1,1,3)
VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)


pdf("VPii.pdf")
par(mfrow=c(1,1),mar=c(5,2,1,1),cex.axis=0.5, cex.lab=0.5)


plotVariancePartitioning(modelsII, VP = VPi, horiz=F, las=2, cex.axis=0.3, cex.names=0.3)
dev.off()

VP

pdf("VP.pdf")
par(mfrow=c(4,2),mar=c(5,2,1,1),cex.axis=0.5, cex.lab=0.5)
for (i in models_f){
VP[[i]] = computeVariancePartitioning(models[[i]],
group = group, groupnames = groupnames)
plotVariancePartitioning(models[[i]], VP = VP[[i]], main=models_fn[i], horiz=F, las=2, cex.axis=0.6, cex.names=0.6)
}
dev.off()

#Current
#VP

#Check colnames
head(modelsII$X)

head(m_pa$X)
head(m_ln_cop$X)

groupnames = c("salinity","habitat")
group = c(1,1,2,2)

VPi= computeVariancePartitioning(modelsII,
group = group, groupnames = groupnames)

VPi_pa= computeVariancePartitioning(m_pa,
group = group, groupnames = groupnames)

VPi_ln_cop= computeVariancePartitioning(m_ln_cop,
group = group, groupnames = groupnames)

VPi<-VPi_ln_cop

#ORDER plots
#run phyloseq taxa filtering (Order)
taxh = subset(tax, (Order %in% tax_to_keep))

taxh2<-taxh[,c("Division","Phylum","Class","Order")]

VPtable<-data.frame(t(VPi$vals))

#Check colnames
colnames(VPtable)<-c("salinity","habitat","substrate","cluster","time","space")

VPtable$Class<-taxh2$Class[match(rownames(VPtable), taxh2$Order)]
VPtable$Phylum<-taxh2$Phylum[match(rownames(VPtable), taxh2$Order)]
VPtable$Division<-taxh2$Division[match(rownames(VPtable), taxh2$Order)]
VPtable$Order<-rownames(VPtable)

#Check colnames
VPtable2<-melt(VPtable, id=c("Division","Phylum","Class","Order"), measure=c("salinity","habitat","time","space","substrate","cluster"))


#Plot variance partition
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc/tax_VP")

for (i in 1:length(unique(VPtable2$Division)))
{
VPtable2.2<-subset(VPtable2, Division==unique(VPtable2$Division)[i])
ggplot(VPtable2.2, aes(Order, value, fill=variable)) + geom_bar(stat="identity") + facet_wrap(Phylum~., scales="free_x") + labs(title="Variance partitioning", x ="Order", y = "variable", fill = "variable") + theme_bw() + theme(axis.text.y = element_text(size=5), axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=5)) + scale_fill_hue(c=45, l=80)
ggsave(paste(as.character(unique(VPtable2$Division)[i]),"barplot_vp_phylum_p25order_scale_log.pdf"))
}


sd_vp2<-ddply(VPtable2, c("Division", "variable"), summarise, 
M = mean(value), SD = sd(value))

ggplot(sd_vp2, aes(x=Division, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Division, scale="free_x")
ggsave("division_vp_p25order.pdf")

sd_vp3<-ddply(VPtable2, c("Division", "Phylum", "variable"), summarise, 
M = mean(value), SD = sd(value))

for (i in 1:length(unique(sd_vp3$Division)))
{
sd_vp3.2<-subset(sd_vp3, Division==unique(sd_vp3$Division)[i])
ggplot(sd_vp3.2, aes(x=Phylum, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Phylum, scale="free_x")
ggsave(paste(as.character(unique(sd_vp3$Division)[i]),"barplot_vp_p25order.pdf"))
}


sd_vp4<-ddply(VPtable2, c("Division","Phylum", "Class", "variable"), summarise, 
M = mean(value), SD = sd(value))

for (i in 1:length(unique(sd_vp3$Division)))
{
sd_vp4.2<-subset(sd_vp4, Division==unique(sd_vp4$Division)[i])
ggplot(sd_vp4.2, aes(x=Class, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Class, scale="free_x")
ggsave(paste(as.character(unique(sd_vp3$Division)[i]),"barplot_vp_p25orderII.pdf"))
}



#CLASS plots
#run phyloseq taxa filtering (Class)

taxh = subset(tax, (Class %in% tax_to_keep))

#If bacteria, fix colnames
colnames(taxh)[1]<-"Division"
#continue..

taxh2<-taxh[,c("Division","Phylum","Class")]

VPtable<-data.frame(t(VPi$vals))

#Import R2
pr2<-cbind(models$spNames,MF$SR2)
colnames(pr2)<-c("Class","SR2")
pr2<-data.frame(pr2)
pr2$SR2<-as.numeric(pr2$SR2)

#Check colnames
colnames(VPtable)<-c("habitat","salinity","substrate","time","space")

VPtable$Phylum<-taxh2$Phylum[match(rownames(VPtable), taxh2$Class)]
VPtable$Division<-taxh2$Division[match(rownames(VPtable), taxh2$Class)]
VPtable$Class<-rownames(VPtable)

VPtable$SR2<-pr2$SR2[match(rownames(VPtable), pr2$Class)]
VPtable$residual<-1-VPtable$SR2

nVPtable<-VPtable[,1:5]*VPtable[,9]
nVPtable$residual<-VPtable$residual
nVPtable[7:9]<-VPtable[,6:8]
colnames(nVPtable)[7:9]<-c("Phylum","Division","Class")

#Check colnames
VPtable2<-melt(nVPtable, id=c("Division","Phylum","Class"), measure=c("habitat","salinity","substrate","time","space","residual"))

VPtable2x<-melt(VPtable, id=c("Division","Phylum","Class"), measure=c("habitat","salinity","substrate","time","space"))



#Plot variance partition
setwd("/home/mpavila/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/both_silva/results/hmsc/tax_VP")


VPtable2.0 <- VPtable2[order(VPtable2$Division, VPtable2$Phylum, VPtable2$Class),]

VPtable2.0$Phylum<-factor(VPtable2.0$Phylum, levels=unique(VPtable2.0$Phylum))
VPtable2.0$Class<-factor(VPtable2.0$Class, levels=unique(VPtable2.0$Class))

VPtable2.0$variable<-factor(VPtable2.0$variable, levels=c("residual","habitat","salinity","substrate","time","space"))

VPtable2.0x <- VPtable2x[order(VPtable2x$Division, VPtable2x$Phylum, VPtable2x$Class),]

VPtable2.0x$Phylum<-factor(VPtable2.0x$Phylum, levels=unique(VPtable2.0x$Phylum))
VPtable2.0x$Class<-factor(VPtable2.0x$Class, levels=unique(VPtable2.0x$Class))

VPtable2.0x$variable<-factor(VPtable2.0x$variable, levels=c("habitat","salinity","substrate","time","space"))

ggplot(data = VPtable2.0, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#999999", "#CC9966", "#006699", "#CC66CC", "#009900", "#666666"),labels = c(
paste("residual",round(mean(nVPtable$residual),digits=2)),
paste("habitat",round(mean(nVPtable$habitat),digits=2)), paste("salinity",round(mean(nVPtable$salinity),digits=2)), paste("substrate",round(mean(nVPtable$substrate),digits=2)), paste("time",round(mean(nVPtable$time),digits=2)), paste("space",round(mean(nVPtable$space),digits=2))))
ggsave("barplot_vp_big_division_p60class_rich_polysal_thin25_s1000_tr02.pdf")


ggplot(data = VPtable2.0x, aes(x = Class, y = value, fill = variable)) + 
geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
facet_grid(~Phylum+Division, switch = "x", scales = "free_x", space = "free_x") +
theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(angle = 90, size = 6), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=6)) + xlab("Clade") +
scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
scale_fill_manual(values=c("#CC9966", "#006699", "#CC66CC", "#009900", "#666666"),labels = c(
paste("habitat",round(mean(VPtable$habitat),digits=2)), paste("salinity",round(mean(VPtable$salinity),digits=2)), paste("substrate",round(mean(VPtable$substrate),digits=2)), paste("time",round(mean(VPtable$time),digits=2)), paste("space",round(mean(VPtable$space),digits=2))))
ggsave("barplot_vp_big_division_p60class_polysal_thin500_s2000_tr02x.pdf")



for (i in 1:length(unique(VPtable2$Division)))
{
VPtable2.2<-subset(VPtable2, Division==unique(VPtable2$Division)[i])
ggplot(VPtable2.2, aes(Class, value, fill=variable)) + geom_bar(stat="identity") + facet_wrap(Phylum~., scales="free_x") + labs(title="Variance partitioning", x ="Order", y = "variable", fill = "variable") + theme_bw() + theme(axis.text.y = element_text(size=5), axis.text.x = element_text(angle = 90, hjust = 1, size=5, vjust=0), strip.text.x = element_text(margin = margin(0.06,0,0.06,0, "cm")), strip.text = element_text(size=5)) + scale_fill_hue(c=45, l=80)
ggsave(paste(as.character(unique(VPtable2$Division)[i]),"barplot_vp_phylum_classp75_thin1.pdf"))
}


sd_vp2<-ddply(VPtable2, c("Division", "variable"), summarise, 
M = mean(value), SD = sd(value))

ggplot(sd_vp2, aes(x=Division, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Division, scale="free_x")
ggsave("division_vp_classp75_root.pdf")

sd_vp3<-ddply(VPtable2, c("Division", "Phylum", "variable"), summarise, 
M = mean(value), SD = sd(value))

for (i in 1:length(unique(sd_vp3$Division)))
{
sd_vp3.2<-subset(sd_vp3, Division==unique(sd_vp3$Division)[i])
ggplot(sd_vp3.2, aes(x=Phylum, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Phylum, scale="free_x")
ggsave(paste(as.character(unique(sd_vp3$Division)[i]),"barplot_vp_phylum_classp75_root.pdf"))
}


sd_vp4<-ddply(VPtable2, c("Division","Phylum", "Class", "variable"), summarise, 
M = mean(value), SD = sd(value))

for (i in 1:length(unique(sd_vp3$Division)))
{
sd_vp4.2<-subset(sd_vp4, Division==unique(sd_vp4$Division)[i])
ggplot(sd_vp4.2, aes(x=Class, y=M, fill=variable)) + geom_bar(stat="identity", colour="black",position="dodge") + theme_bw() + 
geom_errorbar(aes(ymin=M, ymax=M+SD), linewidth=0.5, width=0.25, position=position_dodge(0.9)) + facet_wrap(~Class, scale="free_x")
ggsave(paste(as.character(unique(sd_vp3$Division)[i]),"barplot_vp_classp75_root.pdf"))
}



#Salinity gradient plot
i=104
hab_n<-c("sand","rocks","eelgrass")
pdf("salinity_gradient_test.pdf")
Gradient = constructGradient(models, focalVariable = "Salinity")
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "Y", index = i,showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[1],xlabel = NA,showPosteriorSupport=F)
Gradient = constructGradient(models, focalVariable = "Salinity", 
non.focalVariables = list(habitat = 1))
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "Y", index = i, showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[2],xlabel = NA,showPosteriorSupport=F)
Gradient = constructGradient(models, focalVariable = "Salinity", 
non.focalVariables = list(habitat = 2))
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "Y", index = i, showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[3],xlabel = NA,showPosteriorSupport=F)
dev.off()

i=104
hab_n<-c("sand","rocks","eelgrass")
pdf("salinity_gradient_test_rich.pdf")
Gradient = constructGradient(models, focalVariable = "Salinity")
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "S",showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[1],xlabel = NA,showPosteriorSupport=F)
Gradient = constructGradient(models, focalVariable = "Salinity", 
non.focalVariables = list(habitat = 1))
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "S", showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[2],xlabel = NA,showPosteriorSupport=F)
Gradient = constructGradient(models, focalVariable = "Salinity", 
non.focalVariables = list(habitat = 2))
predY = predict(models, Gradient = Gradient, expected = TRUE)
plotGradient(models, Gradient, pred = predY, measure = "S", showData = TRUE, cex=c(0.08,0.08,0.08), main=hab_n[3],xlabel = NA,showPosteriorSupport=F)
dev.off()



# Run in results folder

library(data.table)
library(dplyr)
library(DECIPHER)
library(seqinr)
library(ape)
library(lattice)

seqs <- fread("COSQ_All_MOTUs_classified.tsv",header=TRUE,sep="\t")

seqtable<-data.frame(sseqid=seqs$sseqid,Seq=seqs$sequence,kingdom=seqs$kingdom,phylum=seqs$phylum,class=seqs$class,order=seqs$order,family=seqs$family,genus=seqs$genus,species=seqs$species)

head(seqtable)

#Remove non-target groups
seqtable<-subset(seqtable,kingdom=="Metazoa")

#Remove duplicate accession numbers. 
seqtable<-seqtable[!duplicated(seqtable$sseqid), ]

seqtable$Seq <- as.character(seqtable$Seq)

#remove gaps
seqtable$Seq <- gsub( "-", "", seqtable$Seq)
#remove duplicates
#seqtable <- seqtable[!duplicated(seqtable), ] #No duplicates present
#remove rows with incomplete taxonomy
seqtable <- seqtable[complete.cases(seqtable),]

#dir.create("barcode_gaps")

write.table(seqtable,"barcode_gaps/COSQ_metazoa_classified_cleaned.txt",col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

#Identify duplicated species
seqtable_duplicated_species <- seqtable[duplicated(seqtable$species), ]
#list duplicated species
duplicated_species_list <- unique(seqtable_duplicated_species$species)

species_identity_table <- c()

for (i in duplicated_species_list){

rows <- seqtable_duplicated_species[which(seqtable_duplicated_species$species == i),]
seqs <- rows$Seq

#if more than 100 rows, sample 100 randomly
if (length(seqs) > 100){ 
seqs <- seqs[sample(length(seqs), 100)] 
}
if (length(seqs) > 1){ 
seqsstring <- DNAStringSet(seqs)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

# calculate identity matrix
y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

min=NA
min75=NA
mean=NA

min <- min(identity_values)
cutoff<- round(length(identity_values)*0.25,0)
if (cutoff == 0){cutoff = 1}
min75 <- identity_values[order(identity_values)][cutoff]
mean <- mean(identity_values)
samples <- nrow(rows)

row <- cbind(rows[1,c(3:8)],min,min75,mean,samples)

species_identity_table <- rbind(species_identity_table,row)
}
}

write.table(species_identity_table,"barcode_gaps/intraspecific.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

histogram <- hist(species_identity_table$min75,breaks=seq(0.6,1,0.01))
t(t(histogram$density))

values <- as.numeric(species_identity_table$min75)
pdf(paste("barcode_gaps/histogram_species_min75.pdf",sep=""),width=12,height=8)
histogram(values, breaks=seq(0.6,1,0.005),col="#99cc00")
dev.off()

###### Genus level (interspecific)

#Get genus list
genus_list <- unique(seqtable$genus)

genus_identity_table <- c()

for (i in genus_list){

rows <- seqtable[which(seqtable$genus == i),]

if (length(unique(rows$species)) > 1){

species_list <- unique(rows$species)

sequence_list <- c()

for (j in species_list){

species_subset <- rows[which(rows$species == j),]
sequence_list <- c(sequence_list,species_subset$Seq)
}

seqsstring <- DNAStringSet(sequence_list)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

min=NA
max=NA
min95=NA
mean=NA

min <- min(identity_values)
max <- max(identity_values)
cutoff<- round(length(identity_values)*0.05,0)
if (cutoff == 0){cutoff = 1}
min95 <- identity_values[order(identity_values)][cutoff]
mean <- mean(identity_values)
samples <- length(species_list)

row <- cbind(rows[1,c(3:7)],max,min,min95,mean,samples)
genus_identity_table <- rbind(genus_identity_table,row)

}
}

histogram <- hist(genus_identity_table$mean,breaks=seq(0.5,1,0.01))

t(t(histogram$density))

write.table(genus_identity_table,"barcode_gaps/interspecific_genus.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

##### Histogram
values <- as.numeric(genus_identity_table$mean)
pdf("barcode_gaps/histogram_genus.pdf",width=12,height=8)
histogram(values, breaks=seq(0.6,1,0.005),col="#99cc00")
dev.off()

###### Family level (between genera)

#Get family list
family_list <- unique(seqtable$family)

family_identity_table <- c()

for (i in family_list){

rows <- seqtable[which(seqtable$family == i),]

if (length(unique(rows$genus)) > 1){

genus_list <- unique(rows$genus)

sequence_list <- c()

for (j in genus_list){

genus_subset <- rows[which(rows$genus == j),]
sequence_list <- c(sequence_list,genus_subset$Seq)
}

seqsstring <- DNAStringSet(sequence_list)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

min=NA
max=NA
min95=NA
mean=NA

min <- min(identity_values)
max <- max(identity_values)
cutoff<- round(length(identity_values)*0.05,0)
if (cutoff == 0){cutoff = 1}
min95 <- identity_values[order(identity_values)][cutoff]
mean <- mean(identity_values)
samples <- length(species_list)

row <- cbind(rows[1,c(3:7)],max,min,min95,mean,samples)
family_identity_table <- rbind(family_identity_table,row)

}
}

histogram <- hist(family_identity_table$mean,breaks=seq(0.4,1,0.01))
t(t(histogram$density))

write.table(family_identity_table,"barcode_gaps/interspecific_family.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

##### Histogram
values <- as.numeric(family_identity_table$mean)
pdf("barcode_gaps/histogram_family.pdf",width=12,height=8)
histogram(values, breaks=seq(0.5,1,0.01),col="#99cc00")

f = hist(family_identity_table$mean,breaks=seq(0.5,1,0.01))
f$density = f$counts/sum(f$counts)*100
g = hist(genus_identity_table$mean,breaks=seq(0.5,1,0.01))
g$density = g$counts/sum(g$counts)*100
h = hist(species_identity_table$min95,breaks=seq(0.5,1,0.01))
h$density = h$counts/sum(h$counts)*100
par(yaxs="i")
hist(rnorm(100)) 

dev.off()

pdf("barcode_gaps/histogram_fgs.pdf",width=12,height=8)
plot(g,col=rgb(1,0,0,0.5),freq=FALSE,ylim=c(0,20),xlab="Average similarity (%)",ylab="Percentage of taxa (%)")
plot(h,col=rgb(0,0,1,0.5),freq=FALSE,add=T)
plot(f,col=rgb(0,1,0,0.5),freq=FALSE,add=T)
dev.off()

###### Order level (between families)

#Get order list
order_list <- unique(seqtable$order)

order_identity_table <- c()

for (i in order_list){

rows <- seqtable[which(seqtable$order == i),]

if (length(unique(rows$family)) > 1){

family_list <- unique(rows$family)

sequence_list <- c()

for (j in family_list){

family_subset <- rows[which(rows$family == j),]
sequence_list <- c(sequence_list,family_subset$Seq)
}

seqsstring <- DNAStringSet(sequence_list)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

min=NA
max=NA
min95=NA
mean=NA

min <- min(identity_values)
max <- max(identity_values)
cutoff<- round(length(identity_values)*0.05,0)
if (cutoff == 0){cutoff = 1}
min95 <- identity_values[order(identity_values)][cutoff]
cutoff2<- round(length(identity_values)*0.95,0)
max95 <- identity_values[order(identity_values)][cutoff2]
mean <- mean(identity_values)
samples <- length(family_list)

row <- cbind(rows[1,c(3:5)],max,min,max95,min95,mean,samples)
order_identity_table <- rbind(order_identity_table,row)

}
}

histogram <- hist(order_identity_table$mean,breaks=seq(0.4,1,0.1))
t(t(histogram$density))

write.table(order_identity_table,"barcode_gaps/interspecific_order.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

##### Histogram
values <- as.numeric(order_identity_table$mean)

pdf("barcode_gaps/histogram_order.pdf",width=12,height=8)
histogram(values, breaks=seq(0.4,1,0.01),col="#99cc00")

e = hist(order_identity_table$mean,breaks=seq(0.4,1,0.01))
e$density = e$counts/sum(e$counts)*100

par(yaxs="i")
hist(rnorm(100)) 
dev.off()
#plot(e,col=rgb(0.5,0,0.5,0.5),freq=FALSE,ylim=c(0,25),xlab="Average similarity (%)",ylab="Percentage of taxa (%)")
#plot(f,col=rgb(0,1,0,0.5),freq=FALSE,add=T)

###### Order level (between families) boxplot

seqtable <- fread("classified_cleaned.txt",header=TRUE,sep="\t")
colnames(seqtable) <- c("GI","Seq","phylum","class","order","family","genus","species")

#Get order list
order_list <- unique(seqtable$order)

boxes <- c()

for (i in order_list){
cat("\nOrder: ",i,"\n")
rows <- seqtable[which(seqtable$order == i),]

if (length(unique(rows$family)) > 1){

family_list <- unique(rows$family)

for (iteration in c(1:50)){

cat(".")

sequence_list <- c()

for (j in family_list){

family_subset <- rows[which(rows$family == j),Seq]
maxiter <- length(family_subset)
random_number <- sample(1:maxiter, 1)
sequence_list <- c(sequence_list,family_subset[random_number])
}

seqsstring <- DNAStringSet(sequence_list)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

box <- cbind(rep(i, length(identity_values)),identity_values)

boxes <- rbind(boxes,box)
}
}
}

boxes <- as.data.frame(boxes)
colnames(boxes) <- c("Order","Similarity")
boxes$Order <- as.character(boxes$Order)
boxes$Similarity <- as.numeric(as.character(boxes$Similarity))

orders <- c("Acanthuriformes","Atheriniformes","Aulopiformes","Beloniformes","Blenniiformes","Carangiformes","Centrarchiformes","Chaetodontiformes","Cichliformes","Clupeiformes","Cypriniformes","Ephippiformes","Gadiformes","Gerreiformes","Gobiiformes","Gonorynchiformes","Istiophoriformes","Kurtiformes","Labriformes","Lutjaniformes","Mugiliformes","Pempheriformes","Perciformes","Pleuronectiformes","Priacanthiformes","Salmoniformes","Scombriformes","Siluriformes","Spariformes","Syngnathiformes","Tetraodontiformes","Uranoscopiformes","Caprimulgiformes","Charadriiformes","Galliformes","Passeriformes","Pelecaniformes","Phoenicopteriformes","Procellariiformes","Sphenisciformes","Carcharhiniformes","Myliobatiformes","Orectolobiformes","Rajiformes","Squaliformes","Cetacea")
boxes2 <- boxes[which(boxes$Order %in% orders),]

pdf("boxplot-orders.pdf",width=8,height=5)
par(mar=c(5,8,4,1)+.1)
boxplot(Similarity~Order,data=boxes2,horizontal=TRUE,las = 1, boxcol="#b1e3e1", medcol="#b1e3e1", whiskcol="#99cc00", staplecol="#99cc00", outpch=20, outcex=0.2, outcol="#cccccc")
##dev.off()

table(boxes2$Order)

```


###### Class level (between orders) boxplot

#Get class list
class_list <- unique(seqtable$class)

boxes <- c()

for (i in class_list){
cat("\nClass: ",i,"\n")

rows <- seqtable[which(seqtable$class == i),]

if (length(unique(rows$order)) > 1){

order_list <- unique(rows$order)

for (iteration in c(1:50)){

cat(".")

sequence_list <- c()

for (j in order_list){

order_subset <- rows[which(rows$order == j),Seq]
maxiter <- length(order_subset)
random_number <- sample(1:maxiter, 1)
sequence_list <- c(sequence_list,order_subset[random_number])
}

seqsstring <- DNAStringSet(sequence_list)
seqs_aligned <- AlignSeqs(seqsstring,verbose = FALSE)
seqs_aligned <- as.character(seqs_aligned)

y <- do.call("rbind", strsplit(seqs_aligned, "")) 
z <- apply(y, 1, function(x) colMeans(x != t(y)) )
identity_matrix <- 1 - z
identity_values <- round(identity_matrix[lower.tri(identity_matrix)],3)

box <- cbind(rep(i, length(identity_values)),identity_values)

boxes <- rbind(boxes,box)
}
}
}

boxes <- as.data.frame(boxes)
colnames(boxes) <- c("Class","Similarity")
boxes$Class <- as.character(boxes$Class)
boxes$Similarity <- as.numeric(as.character(boxes$Similarity))

classes<-c("Actinopteri","Aves","Chondrichthyes","Mammalia")
boxes2 <- boxes[which(boxes$Class %in% classes),]

hist(boxes2$Similarity,breaks=seq(0,1,0.1))
histogram <- hist(boxes2$Similarity,breaks=seq(0,1,0.1))
histogram$density

```

######### Curves

species <- read.table("intraspecific",header=TRUE,sep="\t")
species_values <- species$min95
genus <- read.table("interspecific_genus.txt",header=TRUE,sep="\t")
genus_values <- genus$mean

table <- c()

for (i in seq(0.9,1,0.01)){

species_above_threshold <- length(species_values[which(species_values>=i)])/length(species_values)*100
genus_below_threshold <- length(genus_values[which(genus_values<=i)])/length(genus_values)*100

row <- cbind(i,species_above_threshold,genus_below_threshold)
table <- rbind(table,row)
}




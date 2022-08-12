args = commandArgs(trailingOnly=TRUE)

library(withr)
library(stringr)

tags<-read.table(args[1],header=FALSE,row.names=1)
tags<-tags[!apply(tags == "", 1, all),]   #Not working - not removing final row which is empty
ngsfilter<-data.frame(matrix(NA,nrow=nrow(tags),ncol=6))
colnames(ngsfilter)<-c("#exp","sample","tags","forward_primer","reverse_primer","extra_information")
ngsfilter$"#exp"<-"coastseq_co1"
ngsfilter$no<-seq(1, nrow(ngsfilter), by=1)
ngsfilter$no<-with_options(
c(scipen=999),
str_pad(ngsfilter$no,3,pad="0")
)
ngsfilter$sample<-paste0(args[2],"_sample_",ngsfilter$no)
ngsfilter$tags<-paste0(tags$V2,":",tags$V3)
ngsfilter$forward_primer<-"GGWACWRGWTGRACWNTNTAYCCYCC"
ngsfilter$reverse_primer<-"TANACYTCNGGRTGNCCRAARAAYCA"
ngsfilter$extra_information<-row.names(tags)
ngsfilter<-subset(ngsfilter,select=-no)

###Make metadata file with original names and corresponding mjolnir agnomens
meta<-subset(ngsfilter,select=c(sample,extra_information))
colnames(meta)<-c("mjolnir_agnomens","original_samples")
write.table(meta,file=args[3],sep="\t",row.names=FALSE,quote=FALSE)

ngsfilter<-subset(ngsfilter,select=-extra_information)
write.table(ngsfilter,file=args[4],sep="\t",row.names=FALSE,quote=FALSE)
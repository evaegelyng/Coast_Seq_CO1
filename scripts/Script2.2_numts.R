
library(Biostrings)
library(stringr)

compare.DNA <- function(x,y){
  as.integer(x) == as.integer(y)
}
numts<-function(datas, is_metazoa=FALSE, motu, datas_length)
{
  # compare only mitochondrial genetic code
  mitochondrial_GC <- c(2,3,4,5,7,11,12,14,15,16,17,18)
  # START
  motu_name = motu
  datas$seq<-as.character(datas$seq)
  
  # remove sequences with different length than the seed
  if (sum(datas$id==motu)==0) { # if the seed has been deleted in previous steps take the first more abundant
    motu = datas$id[which(datas$count==max(datas$count,na.rm = TRUE))[1]]
  }
  correct_length <- datas_length[datas$id==motu]
  datas <- datas[datas_length==correct_length,]
  
  # remove misaligned sequences (more than 30 differences between a sequence and 
  # the seed) 
  misaligned_seqs <- c()
  motu_seq <- DNAString(datas$seq[datas$id == motu])
  for (i in 1:dim(datas)[1]) {
    if(sum(!compare.DNA(motu_seq,DNAString(datas$seq[i])))>=30){
      misaligned_seqs <- c(misaligned_seqs,i)
    }
  } 
  if (length(misaligned_seqs)>0) {
    datas <- datas[-misaligned_seqs,]
  }
  
  # look for the best genetic code, this is the one with less stop codons in all sequences
  # the number of codons stop is multiplied by the number of count of the sequence.
  stops<-matrix(NA,dim(datas)[1],20)
  aa_xung<-matrix(NA,dim(datas)[1],20)
  seq<-DNAStringSet(datas$seq)
  seq<-DNAStringSet(seq,start=2,end=nchar(seq[1]))
  
  for (qq in mitochondrial_GC){
    code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[qq]))
    trans<-translate(seq,genetic.code=code)
    
    # for (k in 1:dim(datas)[1]){
    # nstops <- stringr::str_count(as.character(trans[k]),fixed('*'))
    # stops[k,qq] <- nstops * datas$count[k]
    # }
    nstops <- apply(data.frame(as.character(trans)), 1, function(x){stringr::str_count(x,fixed('*'))})
    stops[,qq] <- nstops * datas$count
    
  }
  
  goodcodes<-which(colSums(stops)==min(colSums(stops),na.rm = T))
  
  # if more than one code have been chosen as good code choose the first as the best.
  # However, if the MOTU is a Metazoan and has 313bp length we check the 5 well preserved aa
  # and the code with less changes per read is the chosen. Also remove sequences
  # with changes in thees positions
  
  if (is_metazoa & (correct_length == 313)){  
    
    for (qq in 1:length(goodcodes))
    {
      code<-getGeneticCode(as.character(GENETIC_CODE_TABLE$id[goodcodes[qq]]))
      trans<-translate(seq,genetic.code=code)
      for (k in 1:dim(datas)[1])
      {
        aa<-strsplit(as.character(trans[k]),split="")
        aa<-unlist(aa)
        bad_aa<-0
        if (aa[20]!="H") {bad_aa<-bad_aa+datas$count[k]} # the number of errors counted are the same as the number of counts of the seq.
        if (aa[23]!="G") {bad_aa<-bad_aa+datas$count[k]}
        if (aa[32]!="N") {bad_aa<-bad_aa+datas$count[k]}
        if (aa[81]!="D") {bad_aa<-bad_aa+datas$count[k]}
        if (aa[95]!="G") {bad_aa<-bad_aa+datas$count[k]}
        aa_xung[k,goodcodes[qq]]<-bad_aa
      }
    }
    
    goodcodes <- goodcodes[which(colSums(aa_xung)[goodcodes]==min(colSums(aa_xung)[goodcodes],na.rm = T))]
    bestcode <- goodcodes[1]
    bestcodename <- GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0 | aa_xung[,bestcode]>0
  } else {
    bestcode<-which(colSums(stops)==min(colSums(stops),na.rm = T))[1]
    bestcodename<-GENETIC_CODE_TABLE$name[bestcode]
    goodcodesnames <- GENETIC_CODE_TABLE$name[goodcodes]
    flag <- stops[,bestcode]>0
  }
  
  # numts
  if (sum(flag)>0) {
    numts_seqs <- data.frame('motu' = motu_name, 'id' = datas$id[flag], 
                           'genetic_code' = bestcodename, 
                           'similar_codes' = paste(goodcodesnames, collapse = " | "))
  } else {
    numts_seqs <- c()
  }
  
  # remove numts
  datas <- datas[(flag==FALSE),]
  
  newlist <- list('no_numts_data' = datas, 'numts_seqs' = numts_seqs)
  
  return(newlist)
}

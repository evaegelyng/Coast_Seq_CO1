#!/usr/bin/env Rscript

#script printing out the threshold to filter with obigrep. The output is used as a bash variable and passed to obigrep:
#Example:
# cat tmp/obiuniq_output.fasta | sort -n | uniq -c > tmp/occurrences_of_counts.txt
# threshold=`Rscript --vanilla scripts/threshold_choice.R tmp/occurrences_of_counts.txt` 
# obigrep -p 'count>${threshold}' sample tmp/obiuniq_output.fasta > tmp/obigrep_output.fasta

args = commandArgs(trailingOnly=TRUE)

filename=args[1]

table = read.table(filename, header=FALSE, col.names=c("Occurrences","Value"))

table[,"Variation"] = +Inf

table[1:(dim(table)[1]-1),"Variation"] = table[1:(dim(table)[1]-1),"Occurrences"]/table[2:dim(table)[1],"Occurrences"]

thr = min(which(table["Variation"]<3))
occ = sum(table[(thr+1):(dim(table)[1]),"Occurrences"])
all_occ = sum(table[,"Occurrences"])

#print the threshold, sequences to keep
cat( paste(thr," ",occ," ",all_occ,sep="") )
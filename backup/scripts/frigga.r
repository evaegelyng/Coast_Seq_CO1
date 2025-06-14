args = commandArgs(trailingOnly=TRUE)

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[1])

# FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances

## Script for combining an abundance TSV file from ODIN with a taxonomy-annotated TSV file from THOR.
## The script will have the name of the library (typically 4 characters, e.g. LIBR) as input.
## The taxonomy file must be called LIBR.ecotag.fasta.annotated.tsv
## The abundance file must be called LIBR.SWARM_output.counts.tsv
## The separator characters for both the taxonomy and abundance tsv files must be tabs.
## The abundances file must include sample columns names starting with "sample."
## The output file is a TSV file called LIBR.All_MOTUs.tsv
## FRIGGA deprecates the owi_combine function from previous pipelines (Project Metabarpark, 2016).
## By Owen S. Wangensteen 

mjolnir6_FRIGGA <- function(lib=NULL){
  # sept = separator characters used in taxonomy-annotated file and abundances file, respectively (default: ';;' )
  message("FRYGGA will produce a combined file.")

  infile=paste0(lib,"_classified.tsv")
  abundances=paste0(lib,"_SWARM_output_counts.tsv")
  outfile=paste0(lib,"_All_MOTUs_classified.tsv")

  message("FRYGGA is reading the ecotag-annotated database from THOR...")
  ecotag_db <- read.table(infile,sep="\t",head=T,stringsAsFactors=F)
  message("FRYGGA has read the taxonomy database, with ", nrow(ecotag_db)," total MOTUs.")
  # EES: Replace column name "qseqid" with "id" to allow merging with OTU table (abundance database)
  colnames(ecotag_db)[2]<-"id"

  message("FRYGGA is reading the abundances database from ODIN...")
  abun_db <- read.table(abundances,sep="\t",head=T,stringsAsFactors=F)
  n_samples <- ncol(abun_db[,substr(names(abun_db),1,7)=="sample."])
  message("FRYGGA has read the abundances database, including ", nrow(abun_db)," total MOTUs and ",n_samples," samples.")

  # Merge databases
  db <- merge(ecotag_db,abun_db,by="id")
  #db$sequence <- db$sequence.y
  #db <- db[substr(names(db),1,9)!="sequence."]
  names(db)[substr(names(db),1,7)=="sample."] <- substr(names(db)[substr(names(db),1,7)=="sample."],8,nchar(names(db)[substr(names(db),1,7)=="sample."]))

  # Delete unnecessary columns
  db <- db[,!(names(db) %in% c("definition","ali_length","avg_quality","count","direction","experiment","forward_match","forward_primer","forward_score",
                             "forward_tag","goodali","head_quality","mid_quality","mode","position","reverse_match","reverse_primer","reverse_score","reverse_tag","score","score_norm",
                             "seq_a_deletion","seq_a_insertion","seq_a_mismatch","seq_a_single","seq_ab_match","seq_b_deletion","seq_b_insertion","seq_b_mismatch","seq_b_single",
                             "seq_length_ori","seq_rank","status","tail_quality","rank.1","scientific_name.1","best_identity.1"))]

  # Reorder total_counts and cluster weight columns
  db <- db[,c(1:4,ncol(db)-1,ncol(db)-2,5:ncol(db)-3,ncol(db))]

  write.table(db,outfile,sep="\t",quote=F,row.names=F)
  message("FRYGGA is done. File ", outfile, " written, including ",nrow(db)," MOTUs with ",sum(db$total_reads)," total reads in ",n_samples," samples.")
  message("(",nrow(db[db$total_reads>1,])," non-singletons MOTUs).")
}

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

args = commandArgs(trailingOnly=TRUE)

# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c(args[1])

# Input identifiers for the individual libraries to be used. It should be a 4-character name, matching the information in the ngsfilter files
lib_prefix <- as.character(args[2])

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[2]) # Eva: same as lib_prefixes as each library is run separately in the gwf workflow

####################
# MJOLNIR pipeline #
####################

# Enter number of cores to be used in parallel for RAN and FREYJA. For best performance, the number of libraries to process x the number of cores should be less or equal than the total cores available in the system.
# e.g.: Here we will use 3 cores x 4 libraries = 12 cores. This can be run in any system with 12 cores or more.
cores <- 16  

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames,cores,lib_prefixes=lib_prefix,R1_motif="_1.",R2_motif="_2.") # EES removed "R" in motifs

mjolnir2_FREYJA <- function(lib_prefix="",cores=1,Lmin=299,Lmax=320,lib="",
                            demultiplexed=F,primer_F="GGWACWRGWTGRACWNTNTAYCCYCC",primer_R="TANACYTCNGGRTGNCCRAARAAYCA",R1_motif="_R1",R2_motif="_R2",obipath=""){
  message("FREYJA will do paired-end alignment, demultiplexing and length filter.")  
  suppressPackageStartupMessages(library(parallel))
  no_cores <- cores*length(lib_prefix)
  old_path <- Sys.getenv("PATH")
  clust <- makeCluster(no_cores)
  X <- NULL
  libslist <- NULL

  if (!demultiplexed){
    for (i in 1:cores) for (j in 1:length(lib_prefix)) {
      #if (cores<10) {formatted_i <- i} else {formatted_i <- sprintf("%02d",i)}
      formatted_i <- sprintf("%02d",i)
      X <- c(X,paste0("illuminapairedend -r ",lib_prefix[j],"_R2_part_",formatted_i,".fastq ",lib_prefix[j],"_R1_part_",formatted_i,".fastq | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",lib_prefix[j],".tsv | obigrep -p \'seq_length>",Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\' --fasta-output > ",lib_prefix[j],".filtered_length_part_",sprintf("%02d",i),".fasta"))
    }  
    for (prefix in lib_prefix) libslist <- paste0(libslist,prefix,".filtered_length_part*.fasta ")
  } else {
      metadata <- read.table(paste0(lib,"_metadata.tsv"),sep="\t",header=T)
      fastqR1_list <- metadata$fastq_name_R1
      agnomens <-  metadata$mjolnir_agnomens
      # Create ngsfilter files
      for (ag in agnomens) writeLines(paste(lib,ag,":",primer_F,primer_R,sep="\t"),paste0("ngsfilter_",ag,".tsv"))
      # Create obitool commands
      for (i in 1:length(agnomens)) {
         X <- c(X,paste0("illuminapairedend -r ",gsub(R1_motif,R2_motif,fastqR1_list[i]), " ",fastqR1_list[i]," | obigrep -p \'score>40.00\' | ngsfilter -t ngsfilter_",agnomens[i],".tsv | obigrep -p \'seq_length>",Lmin,"\' -p \'seq_length<",Lmax,"\' -s \'^[ACGT]+$\' -p \'forward_tag!=None\' -p \'reverse_tag!=None\' --fasta-output > ",agnomens[i],".fasta"))
         libslist <- paste0(libslist,agnomens[i],".filtered_length_part_",sprintf("%02d",i),".fasta ")
      }
  }
    clusterExport(clust, list("X","old_path","obipath"),envir = environment())
    clusterEvalQ(clust, {Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))}) 
    parLapply(clust,X, function(x) system(x,intern=T,wait=T))
    stopCluster(clust)
    # If not demultiplexed, then join all parts into a joined file and then split it into samples
    if (!demultiplexed){
      message("FREYJA is joining filtered reads into a single file.")
      system(paste0("cat ",libslist," > ",lib,".joined.fasta"),intern=T,wait=T)
      message("File ",lib,".joined.fasta written.")
      message("FREYJA will create individual files for each sample.")
      system(paste0("obisplit -t sample ",lib,".joined.fasta"),intern=T,wait=T)
    }
    message("FREYJA is done.")
}

# FREYJA will do the paired-end alignment, demultiplexing & length filtering. It will give individual filtered sample files as an output.
mjolnir2_FREYJA(lib_prefix,cores,Lmin=299,Lmax=320)

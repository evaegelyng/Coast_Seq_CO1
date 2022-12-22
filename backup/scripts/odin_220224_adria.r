args = commandArgs(trailingOnly=TRUE)

# Input path to the OTU table for the final combined library
lib_in <- as.character(args[1])

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[2])

cores <- 32
 
# ODIN: OTU Delimitation Inferred by Networks

## ODIN performs MOTU clustering and (optionally) entropy denoising and it is one of the main steps of MJOLNIR.
## Clustering is performed using SWARM, which will produce MOTUs as networks of unique sequences 
## After SWARM, ODIN will recalculate the abundances for every MOTU in every sample.
## Then (optionally) ODIN proceeds with within-MOTU denoising, using the DnoisE entropy-ratio algorithm for coding regions to get an ESV table. 
## Two obligatory arguments are needed: the name of the library, typically 4 characters, and the number of computing cores.
## Three optional parameters: the clustering distance d (default=13), 
## min_reads_MOTU is the minimum number of reads to keep a MOTU in the MOTU table (default=2),
## and the minimum number of reads to keep an ESV in the final ESV file (default=2).
## Two boolean parameters can be selected: run_swarm can be set as FALSE to save time if a SWARM output is already available.
## And generate_ESV=TRUE (default) will use the DnoisE algorithm to produce an ESV table, along with the MOTU table.
## ODIN deprecates the previous owi_recount_swarm script used in old metabarcoding pipelines (e.g. Project Metabarpark 2015).
## By Owen S. Wangensteen

## Eva Sigsgaard: Change the path to dnoise executable and the dnoise options (remove/correct old (?) options, and add modal length of 313 bp and parallel computing). Move singleton removal to separate workflow target. 

mjolnir4_ODIN_eva <- 
function(lib,cores,d=13,min_reads_MOTU=2,min_reads_ESV=2,run_swarm=TRUE,generate_ESV=FALSE,obipath=""){

  dnoise_path <- "/home/evaes/miniconda3/pkgs/dnoise-1.0-py38_0/lib/python3.8/site-packages/src/DnoisE.py"    # Change this to where the Dnoise executable is
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, obipath, sep = ":"))
  if (run_swarm){
    message("ODIN will cluster sequences into MOTUs with SWARM.")
    system(paste0("swarm -d ",d," -z -t ",cores," -o ",lib,"_SWARM_output -s ",lib,"_SWARM",d,"nc_stats -w ",lib,"_SWARM_seeds.fasta ",lib_in,".vsearch.fasta"),intern=T,wait=T)
    message("ODIN will recount abundances for every MOTU after Swarm.")
  }
  fileswarm=paste0(lib,"_SWARM_output")
  filetab=paste0(lib_in,"_new.tab")
  outfile_MOTU <-paste(fileswarm,"_counts.tsv",sep="")
  outfile_ESV <-paste(fileswarm,"_ESV.tsv",sep="")
  
  get_swarm_size <- function(cadena="="){
    # it gets the positon of the '=' and gets the items after it
    return(as.numeric(gsub(";","",substr(cadena,gregexpr("=",cadena)[[1]][[1]]+1,nchar(cadena)))))
  }
  
  message("ODIN will recount abundances for every MOTU after Swarm.")
  
  # Read cluster list database
  message("1. ODIN is reading SWARM results...")
  swarm_db <- readLines(fileswarm)
  total_swarms <- length(swarm_db)
  message("2. ODIN has read ", total_swarms," total MOTUs.")
  
  message(paste("3. ODIN will now perform a remove of MOTUs with sequences of less than",
                min_reads_MOTU,"reads."))
  # Calculate reads in each cluster and reduce the dataset if min_reads_MOTU>0
  clusters <- strsplit(swarm_db,"; ")
  if (min_reads_MOTU>0) {
    # first reduction of a proportion of MOTUs
    if (min_reads_MOTU>9) i <- 9 else i <- min_reads_MOTU-1
    clusters <- clusters[!(grepl(paste0('size=[0-',i,']'),clusters) & lengths(clusters)==1)]
    cluster_reads <- NULL
    # second reduction
    cluster_reads <- mclapply(clusters,function(x) sum(as.numeric(lapply(X=(x),FUN=get_swarm_size))), mc.cores = cores)
    # for (i in 1:length(clusters)) cluster_reads[i] <- sum(as.numeric(lapply(X=(clusters[[i]]),FUN=get_swarm_size)))
    clusters <- clusters[cluster_reads>=min_reads_MOTU]
    total_swarms <- length(clusters)
  }
  
  message("4. ODIN will now calculate the number of reads in every sample for each MOTU.")
  clusters <- mclapply(clusters,function(x){sub(";.*","",x)}, mc.cores = cores)
  names(clusters) <- mclapply(clusters, function(x) x[[1]], mc.cores = cores)
  
  message("5. ODIN kept only ", total_swarms," MOTUs of size greater than or equal to ",min_reads_MOTU," reads.")
  motu_seqs_names <- stack(clusters) %>% rename(id = values, MOTU = ind)
  # motu_seqs_names <- unlist(clusters, use.names=F)
  
  # # Generate a file with the list of ids of non-singleton clusters
  # motulist <- file(paste0(lib,"_non_singleton_motu_list.txt"),"wt")
  # writeLines(id,motulist)
  # message("ODIN has created the file ",paste0(lib,"_non_singleton_motu_list.txt")," with the list of identifiers of non-singleton MOTUs.")
  
  # Read counts database and keep only the needed clusters
  message("6. ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
  db <- read.table(filetab,sep="\t",head=T)
  numseqs <- nrow(db)
  db <- merge(motu_seqs_names2,db,by="id")
  # db <- db[db$id %in% motu_seqs_names,]
  numseqs_reduced <- nrow(db)
  samples <- sum(grepl('sample',names(db)))
  message("7. ODIN finished reading the Database, which includes ", numseqs," total unique sequences and ",samples," samples.\n",
          "ODIN kept only ", numseqs_reduced," sequences for calculations.")
  
  message("4. ODIN will now calculate the number of reads in every sample for each MOTU.")
  db.total <- split(db[,grepl('sample',names(db))], db$MOTU)
  db.total <- mclapply(db.total,function(x)as.data.frame(t(as.matrix(c(count=sum(x),colSums(x),CLUST_WEIGHT=dim(x)[1])))), mc.cores = cores)
  db.total <- do.call(rbind,db.total)
  db.total <- cbind(data.frame(id=rownames(db.total)), db.total)
  db.total <- merge(db.total,db[,grepl('id|sequence',names(db))],by = "id")
  
  # order the columns
  col_order <- c('id', 'count', 'MOTU', names(db)[grepl('sample',names(db))], 'sequence' )
  db <- db[,col_order]
  s_opt <- min(grep('sample',names(db)))
  z_opt <- max(grep('sample',names(db)))
  
  # print datasets
  write.table(db.total,outfile_MOTU,sep="\t",quote=F,row.names=F)
  write.table(db,outfile_ESV,sep="\t",quote=F,row.names=F)
  
  # divide dataset into different files for THOR
  db.total <- db.total[,c('id','sequence')]
  db.total <- paste(paste0('>',db.total$id),db.total$sequence,sep='\n')
  writeLines(paste0(db.total[[part]],collapse = '\n'),paste0(lib,"_more_than_",min_reads_MOTU,".fasta"))
  
  
  message("File ", lib,"_more_than_",min_reads_MOTU,".fasta", " written")
  
  
}

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN_eva(lib,cores,d=13,run_swarm=TRUE,generate_ESV=FALSE)

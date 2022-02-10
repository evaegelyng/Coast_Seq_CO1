args = commandArgs(trailingOnly=TRUE)

# Input path to the OTU table for the final combined library
lib_in <- as.character(args[1])

# Input name for the final combined library (should be a 4-character name)
lib <- as.character(args[2])

cores <- 32
  
# Adapt ODIN function. Change the path to dnoise executable and the dnoise options (remove/correct old (?) options, and add modal length of 313 bp and parallel computing). Move singleton removal to separate workflow target. 
mjolnir4_ODIN_eva <- function(lib, cores, d = 13, min_reads_MOTU = 2, min_reads_ESV = 2,
    run_swarm = TRUE, generate_ESV = FALSE)
{
    dnoise_path <- "/home/evaes/miniconda3/pkgs/dnoise-1.0-py38_0/lib/python3.8/site-packages/src/DnoisE.py"    # EES changed path to dnoise
    if (run_swarm) {
        message("ODIN will cluster sequences into MOTUs with SWARM.")
        system(paste0("swarm -d ", d, " -z -t ", cores, " -o ",
            lib, ".SWARM_output -s ", lib, ".SWARM", d, "nc_stats -w ",
            lib, ".SWARM_seeds.fasta ", lib_in, ".vsearch.fasta"),
            intern = T, wait = T)
        message("ODIN will recount abundances for every MOTU after Swarm.")
    }
    fileswarm = paste0(lib, ".SWARM_output")
    filetab = paste0(lib_in, ".new.tab")
    outfile <- paste(fileswarm, ".counts.csv", sep = "")
    outfile_ESV <- paste(fileswarm, ".ESV.csv", sep = "")
    get_swarm_size <- function(cadena = "=") {
        return(as.numeric(substr(cadena, gregexpr("=", cadena)[[1]][[1]] +
            1, nchar(cadena))))
    }
    message("ODIN is reading SWARM results...")
    swarm_db <- readLines(fileswarm)
    total_swarms <- length(swarm_db)
    message("ODIN has read ", total_swarms, " total MOTUs.")
    message("ODIN will now calculate the number of reads in every sample for each MOTU.")
    clusters <- strsplit(swarm_db, "; ")
    for (i in 1:total_swarms) for (j in 1:length(clusters[[i]])) if (substr(clusters[[i]][[j]],
        nchar(clusters[[i]][[j]]), nchar(clusters[[i]][[j]])) ==
        ";") {
        clusters[[i]][[j]] <- substr(clusters[[i]][[j]], 1, nchar(clusters[[i]][[j]]) -
            1)
    }
    cluster_reads <- NULL
    for (i in 1:total_swarms) cluster_reads[i] <- sum(as.numeric(lapply(X = (clusters[[i]]),
        FUN = get_swarm_size)))
    swarm_db_reduced <- swarm_db[cluster_reads >= min_reads_MOTU]
    clusters <- strsplit(swarm_db_reduced, "; ")
    total_swarms_reduced <- length(swarm_db_reduced)
    id <- NULL
    for (i in 1:total_swarms_reduced) for (j in 1:length(clusters[[i]])) {
        clusters[[i]][[j]] <- sub(";.*", "", clusters[[i]][[j]])
        id[i] <- clusters[[i]][1]
    }
    names(clusters) <- id
    message("ODIN kept only ", total_swarms_reduced, " MOTUs of size greater than or equal to ",
        min_reads_MOTU, " reads.")
    necesarios <- unlist(clusters, use.names = F)
    motulist <- file(paste0(lib, "_non_singleton_motu_list.txt"),
        "wt")
    writeLines(id, motulist)
    message("ODIN has created the file ", paste0(lib, "_non_singleton_motu_list.txt"),
        " with the list of identifiers of non-singleton MOTUs.")
    message("ODIN is reading the abundance database. This could take Him a while, since He has just one eye left, after all.")
    db <- read.table(filetab, sep = "\t", head = T)
    numseqs <- nrow(db)
    db$id <- gsub(";", "", db$id)
    db <- db[db$id %in% necesarios, ]
    numseqs_reduced <- nrow(db)
    samples <- length(names(db)[substr(names(db), 1, 6) == "sample"])
    message("ODIN finished reading the Database, which includes ",
        numseqs, " total unique sequences and ", samples, " samples.")
    message("ODIN kept only ", numseqs_reduced, " sequences for calculations.")
    db.total <- merge(data.frame(id), db, by = "id")
    id <- db.total$id
    numclust <- nrow(db.total)
    if (generate_ESV)
        dir.create("MOTU_tsv", showWarnings = FALSE)
    for (fila in 1:numclust) {
        head <- id[fila]
        tails <- unlist(clusters[names(clusters) == head])
        db.reduced <- db[db$id %in% tails, ]
        if (generate_ESV)
            write.table(db.reduced, paste0("MOTU_tsv/", head),
                sep = "\t", quote = F, row.names = F)
        suma <- colSums(db.reduced[, substr(names(db.total),
            1, 6) == "sample"])
        db.total[fila, substr(names(db.total), 1, 6) == "sample"] <- suma
        db.total$cluster_weight[fila] <- nrow(db.reduced)
        message("Cluster ", fila, " / ", numclust, " ready, including ",
            db.total$cluster_weight[fila], " sequences.", "\r",
            appendLF = FALSE)
    }
    db.total$total_reads <- rowSums(db.total[, substr(names(db.total),
        1, 6) == "sample"])
    names(db.total[substr(names(db.total), 1, 6) == "sample"]) <- substr(names(db.total[substr(names(db.total),
        1, 6) == "sample"]), 8, nchar(names(db.total[substr(names(db.total),
        1, 6) == "sample"])))
    write.table(db.total[, c(1:(ncol(db.total) - 3), (ncol(db.total) -
        1):ncol(db.total), (ncol(db.total) - 2))], outfile, sep = ";",
        quote = F, row.names = F)
    message("File ", outfile, " written")
    if (generate_ESV) {
        message("ODIN will generate now a list of ESVs for every non-singleton MOTU, using DnoisE.")
        sample_cols <- (1:ncol(db.total))[substr(names(db.total),
            1, 6) == "sample"]
        start_samp <- sample_cols[1]
        end_samp <- sample_cols[length(sample_cols)]
        suppressPackageStartupMessages(library(parallel))
        clust <- makeCluster(cores)
        X <- NULL
        for (motu in id) X <- c(X, paste0("python3 ", dnoise_path, " --csv_input MOTU_tsv/",
            motu, " --csv_output ", motu, " -s ", start_samp,
            " -z ", end_samp, " -n 'count' -p 1 -y T -c 16 -m 313"))
        clusterExport(clust, "X", envir = environment())
        parLapply(clust, X, function(x) system(x, intern = T,
            wait = T))
        stopCluster(clust)
        message("ODIN will now merge all ESVs into a final ESV table.")
        ESV_tables <- NULL
        for (i in 1:length(id)) {
            ESV_tables[[i]] <- read.table(paste0(id[i], "_Adcorr_denoised_ratio_d.csv"),
                sep = ",", head = T)
            ESV_tables[[i]]$definition <- id[i]
            ESV_tables[[i]] <- ESV_tables[[i]][ESV_tables[[i]]$count >=
                min_reads_ESV, ]
        }
        ESV_table <- ESV_tables[[1]]
        for (i in 2:length(id)) ESV_table <- rbind(ESV_table,
            ESV_tables[[i]])
        write.table(ESV_table, outfile_ESV, sep = ";", quote = F,
            row.names = F)
        message("File ", outfile_ESV, " written with ", nrow(ESV_table),
            " ESVs in ", length(id), " MOTUs")
    }
        message("ODIN is now removing singleton MOTUs from the fasta output file, to make THOR's work an easier task.")
        system(paste0("sed -i 's/;size/ size/g' ", lib, ".SWARM_seeds.fasta"),
        intern = T, wait = T)
}

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN_eva(lib,cores,d=13,generate_ESV=TRUE)

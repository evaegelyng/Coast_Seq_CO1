args = commandArgs(trailingOnly=TRUE)

fila = args[1]

db<-read.table("db.tsv",sep="\t",quote=F,row.names=F)
clusters <- readRDS("clusters.rds")
id <- readRDS("id.rds")
db.total <- merge(data.frame(id),db,by="id") # This will keep just the heads

head <- id[fila]
tails <- unlist(clusters[names(clusters)==head])
db.reduced <- db[db$id %in% tails,]
suma <- colSums(db.reduced[,substr(names(db.total),1,6)=="sample"])
db.total[fila,substr(names(db.total),1,6)=="sample"] <- suma
db.total$cluster_weight[fila] <- nrow(db.reduced)
db.fila<-db.total[fila,]
write(db.fila)
message("Cluster ", fila, " / ",numclust, " ready, including ", db.total$cluster_weight[fila]," sequences.","\r",appendLF = FALSE)

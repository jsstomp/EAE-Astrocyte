file_list <- list.files("DEA_Results", full.names=T)
all_results <- lapply(file_list, read.delim, sep= "\t", header=T, row.names=1)
all_results <- do.call(rbind, all_results)
tt <- head(all_results[order(all_results$padj),], 10)
write.table(tt, "Results/top_table.txt",quote=F,col.names=NA,row.names=T,sep="\t")
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
IDs <- rownames(all_results)
G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol", "description"), values = IDs, mart = mart)
all_results$ensembl_gene_id <- IDs
m <- merge(all_results, G_list, by="ensembl_gene_id")
m <- m[-9]
dim(m)
head(m[order(m$padj),], 10)

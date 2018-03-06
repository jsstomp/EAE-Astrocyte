file_list <- list.files("DEA_Results", full.names=T)
options(stringsAsFactors = F)
all_results <- lapply(file_list, read.delim, sep= "\t", header=T)
all_results <- do.call(rbind, all_results)
aa <- all_results[order(all_results$ensembl_gene_id, abs(all_results$padj) ), ] #sort by id and reverse of abs(value)
all_results <- aa[ !duplicated(aa$ensembl_gene_id), ]              # take the first row within each id



# To make sure that the gene symbols of the Ensembl id's are known,
# a txt file is used of a gtf file (used during the alignment).
# geneExp <- read.csv(file = "/Users/marissadubbelaar/Desktop/Genomes/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)
geneExp <- read.csv(file = "RawFiles/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)

# Only the gene id (ensemble id) and the gene name (genesymbol) was
# obtained from this file.
geneNames <- cbind(gsub("gene_id ", "", as.matrix(geneExp[,1])), gsub("gene_name ", "", as.matrix(geneExp[,5])), trimws(apply(geneExp, 1, function(x) tail(na.omit(x), 1))))

gene_names <- as.data.frame(unique(geneNames[,1:2]))
colnames(gene_names) <- c("ensembl_gene_id", "external_gene_name")
m <- merge(all_results, gene_names, by="ensembl_gene_id")
head(m2)
head(all_results)
library(pheatmap)
library(edgeR)
library(RColorBrewer)
mycols <- c('#1f78b4','#e31a1c','#33a02c','#696969','#b2df8a','#fdbf6f',
             '#ff7f00','#b15928','#a6cee3','#cab2d6','#fb9a99','#6a3d9a')
countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
countData <- countData[which(rownames(countData) %in% m$ensembl_gene_id), ]
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,c(4:6,13)]
colData$RP <- paste(colData$Region, colData$Population, sep = "_")
e <- cpm(countData)
cols <- palette(mycols)[as.fumeric(colData$groupID)]
annotation <- data.frame(Group = colData$groupID, Region = colData$RP)
rownames(annotation) <- colnames(e)
mycols <- list(Group = c(C_HB_A = 'darkviolet', E4_HB_A = 'darkorchid', Ech_HB_A = 'deeppink', E1_HB_A = 'hotpink', C_HB_AG = 'navy', E1_HB_AG = 'royalblue',
                         E4_HB_AG = 'darkcyan', Ech_HB_AG = 'deepskyblue', C_SC_A = 'yellow', E1_SC_A = 'gold', E4_SC_A = 'gold3', Ech_SC_A = 'orangered'),
               cluster = (c("cadetblue" , "hotpink", "navy", "orangered", "gold", "darkviolet")))
hmcol <- colorRampPalette(c("magenta", "dodgerblue2", "white", "tan3","firebrick3"))(n=99)


heatmap.2(e, scale = "row", trace = "none", ColSideColors = cols, col = hmcol, cexRow = 0.05, hclustfun = "euclidian" )
pheatmap(e, scale = "row", trace = "none", annotation_col = as.character(colData$groupID), annotation_colors = mycols[1],annotation_legend = FALSE, col = hmcol, fontsize_row = 0.05, fontsize_col = 6)
?pheatmap

pres <- pheatmap(e, scale = "row", annotation = annotation, annotation_colors = mycols, col = hmcol, fontsize_row = 0.05, fontsize_col = 6)
pres.clust <- cbind(pres, cluster = cutree(pres$tree_row, k=6))
pres.clust <- as.data.frame(unlist(pres.clust[,2]))
colnames(pres.clust) <- c("cluster")
ann.cols <- list(cluster=colorRampPalette(c("dodgerblue2", "firebrick3")))

pheatmap(e, scale = "row", show_rownames = F, show_colnames = F, cutree_rows = 6, annotation = annotation,
         annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)
e2 <- e[,order(colData$RP, colData$group)]
pheatmap(e2, scale = "row", show_rownames = F, show_colnames = F, cutree_rows = 6, cluster_cols=FALSE,
         annotation = annotation, annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)

pres.clust$ensembl_gene_id <- rownames(pres.clust)
m3 <- merge(m, pres.clust, by = "ensembl_gene_id")
m3 <- m3[order(m3$cluster),]
rownames(m3) <- m3$external_gene_name
m3 <- m3[c(1:7,9)]

tt <- head(m[order(m$padj),], 10)
write.table(tt, "Results/top_table.txt",quote=F,col.names=NA,row.names=T,sep="\t")
write.csv(m3, file="cluster_metascape_input.csv")

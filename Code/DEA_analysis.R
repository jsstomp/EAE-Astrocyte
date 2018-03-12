file_list <- list.files("DEA_Results", full.names=T)
HBAG_list <- grep("HBAG_.+HBAG\\.", file_list, value=TRUE)
HBA_list <- grep("HBA_.+HBA\\.", file_list, value=TRUE)
SCA_list <- grep("SCA_.+SCA\\.", file_list, value=TRUE)
options(stringsAsFactors = F)

library(pheatmap)
library(edgeR)
library(RColorBrewer)




TBN <- function(filenames, region) {
  # load and bind given files
  result <- lapply(filenames, read.delim, sep="\t", header=T)
  result <- do.call(rbind, result)
  # sort by adjusted p-value ascending and remove duplicates with least significant adjusted p-value
  aa <- result[order(result$ensembl_gene_id, abs(result$padj) ), ]  #sort by id ascending from absolute of adjusted p-value
  result <- aa[ !duplicated(aa$ensembl_gene_id), ]  #take first row within each id
  # merge result with genename dataset to get the corrosponding gene name for each ensembl id
  m <- merge(result, gene_names, by="ensembl_gene_id")
  # perpare countdata for heatmapping
  if(region=="all"){
    count_data <- countData
  }
  else {
    count_data <- countData[,grepl(region, colnames(countData))]
  }
  count_data <- count_data[which(rownames(count_data) %in% m$ensembl_gene_id), ]
  e <- cpm(count_data)
  # prepare coldata
  col_data <- colData[which(rownames(colData) %in% colnames(count_data)),]
  #heatmap
  annotation <- data.frame(Group = col_data$groupID, Region = col_data$RP)
  rownames(annotation) <- colnames(e)
  pres <- pheatmap(e, scale = "row", annotation = annotation, annotation_colors = mycols, col = hmcol, fontsize_row = 0.05, fontsize_col = 6)
  pres.clust <- cbind(pres, cluster = cutree(pres$tree_row, k=6))
  pres.clust <- as.data.frame(unlist(pres.clust[,2]))
  colnames(pres.clust) <- c("cluster")
  ann.cols <- list(cluster=colorRampPalette(c("dodgerblue2", "firebrick3")))
  print(colnames(pres.clust))
  #rownames(annotation) <- colnames(e)
  pheatmap(e, scale = "row", show_rownames = F, show_colnames = F, cutree_rows = 6, annotation = annotation,
           annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)
  e2 <- e[,order(col_data$RP, col_data$group)]
  pheatmap(e2, scale = "row", show_rownames = F, show_colnames = F, cluster_cols=FALSE, annotation = annotation,
           annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)
  
}

# To make sure that the gene symbols of the Ensembl id's are known,
# a txt file is used of a gtf file (used during the alignment).
# geneExp <- read.csv(file = "/Users/marissadubbelaar/Desktop/Genomes/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)
geneExp <- read.csv(file = "RawFiles/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)), fill = T, na.strings=c("","NA"), stringsAsFactors = F)

# Only the gene id (ensemble id) and the gene name (genesymbol) was
# obtained from this file.
geneNames <- cbind(gsub("gene_id ", "", as.matrix(geneExp[,1])), gsub("gene_name ", "", as.matrix(geneExp[,5])), trimws(apply(geneExp, 1, function(x) tail(na.omit(x), 1))))

gene_names <- as.data.frame(unique(geneNames[,1:2]))
colnames(gene_names) <- c("ensembl_gene_id", "external_gene_name")

mycols <- list(Group = c(C_HB_A = 'darkviolet', E4_HB_A = 'darkorchid', Ech_HB_A = 'deeppink', E1_HB_A = 'hotpink', C_HB_AG = 'navy', E1_HB_AG = 'royalblue',
                         E4_HB_AG = 'darkcyan', Ech_HB_AG = 'deepskyblue', C_SC_A = 'yellow', E1_SC_A = 'gold', E4_SC_A = 'gold3', Ech_SC_A = 'orangered'),
               cluster = (c("cadetblue" , "hotpink", "navy", "orangered", "gold", "darkviolet")))
hmcol <- colorRampPalette(c("magenta", "dodgerblue2", "white", "tan3","firebrick3"))(n=99)

countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,c(4:6,13)]
colData$RP <- paste(colData$Region, colData$Population, sep = "_")


count <- 0
for(files in list(file_list,HBA_list,HBAG_list,SCA_list)){
  if(count==1){
    region = "HB_A_"
  }
  else if(count==2){
    region = "HB_AG_"
  }
  else if(count==3){
    region = "SC_A_"
  }
  else{
    region = "all"
  }
  count = count + 1
  TBN(files, region )
}
# 
# all_results <- lapply(file_list, read.delim, sep= "\t", header=T)
# all_results <- do.call(rbind, all_results)
# aa <- all_results[order(all_results$ensembl_gene_id, abs(all_results$padj) ), ] #sort by id and reverse of abs(value)
# all_results <- aa[ !duplicated(aa$ensembl_gene_id), ]              # take the first row within each id
# 
# m <- merge(all_results, gene_names, by="ensembl_gene_id")
# head(m2)
# head(all_results)
# 
# 
# cols <- palette(mycols)[as.fumeric(colData$groupID)]
# annotation <- data.frame(Group = colData$groupID, Region = colData$RP)
# rownames(annotation) <- colnames(e)
# 
# heatmap.2(e, scale = "row", trace = "none", ColSideColors = cols, col = hmcol, cexRow = 0.05, hclustfun = "euclidian" )
# pheatmap(e, scale = "row", trace = "none", annotation_col = as.character(colData$groupID), annotation_colors = mycols[1],annotation_legend = FALSE, col = hmcol, fontsize_row = 0.05, fontsize_col = 6)
# ?pheatmap
# 
# 
# pres.clust$ensembl_gene_id <- rownames(pres.clust)
# m3 <- merge(m, pres.clust, by = "ensembl_gene_id")
# m3 <- m3[order(m3$cluster),]
# rownames(m3) <- m3$external_gene_name
# m3 <- m3[c(1:7,9)]

tt <- head(m[order(m$padj),], 10)
write.table(tt, "Results/top_table.txt",quote=F,col.names=NA,row.names=T,sep="\t")
write.csv(m3, file="cluster_metascape_input.csv")







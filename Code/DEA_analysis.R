####################################################################
# Author: Jafta Stomp
# Date: 15-02-2018)
# Description: 
#   This script analyzes results of the differential expression
#   anaylis run in DESeq2, also makes heatmaps of different regions.
####################################################################

####################################################################
#                           IMPORTS                                #
####################################################################
file_list <- list.files(paste("Results/", prefix,"/DEA_Results/", sep=""), full.names=T)
HBAG_list <- grep("HBAG_.+HBAG\\.", file_list, value=TRUE)
HBA_list <- grep("HBA_.+HBA\\.", file_list, value=TRUE)
SCA_list <- grep("SCA_.+SCA\\.", file_list, value=TRUE)
options(stringsAsFactors = F)

cat("Loading required libraries.\n")
suppressMessages(library(pheatmap))
suppressMessages(library(edgeR))
suppressMessages(library(RColorBrewer))
suppressMessages(library(fpc))

####################################################################
#                           FUNCTIONS                              #
####################################################################
cat("Building functions.\n")
analyzer <- function(filenames, region) {
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
  # heatmap color annotation
  annotation <- data.frame(Group = col_data$groupID, Region = col_data$RP)
  rownames(annotation) <- colnames(e)
  # calculate the number of clusters using the Calinski-Harabasz index
  nb <- pamk(e, krange = 3:8, criterion = "ch", usepam = T)
  nc <- nb$nc
  # make an initial pheatmap object to use for defining clusters of genes
  pdf(file = NULL) # null pdf to avoid pheatmap from sending this to a file
  pres <- pheatmap(e, scale = "row", annotation = annotation, annotation_colors = mycols, col = hmcol, fontsize_row = 0.05, fontsize_col = 6)
  invisible(dev.off())
  # cut the pheatmap tree using cutree and the predetermined number of clusters
  pres.clust <- cbind(pres, cluster = cutree(pres$tree_row, k=nc))
  pres.clust <- as.data.frame(unlist(pres.clust[,2]))
  # Set clusters to colorcodes as requested by biologist
  result$cluster <- pres.clust[order(rownames(pres.clust)),]
  result$cluster[result$cluster==1] <- "#ff5300"
  result$cluster[result$cluster==2] <- "#3f5eba"
  result$cluster[result$cluster==3] <- "#e9003a"
  result$cluster[result$cluster==4] <- "#56c7e9"
  result$cluster[result$cluster==5] <- "#006633"
  result$cluster[result$cluster==6] <- "#955659"
  result$cluster[result$cluster==7] <- "#14e337"
  result$cluster[result$cluster==8] <- "#ff73f9"
  # Write results to csv file (per region) to be used by biologist
  write.csv(result, file = paste("Results/",prefix,"/de_genes_",gsub("_","",region), sep=""))
  colnames(pres.clust) <- c("cluster")
  ann.cols <- list(cluster=colorRampPalette(c("dodgerblue2", "firebrick3")))
  # make unsupervised heatmap
  pheatmap(e, main = paste("Heatmap Unsupervised Clustering", gsub("_","",region)), scale = "row", show_rownames = F, show_colnames = F, cutree_rows = nc, annotation = annotation,
           annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)
  e2 <- e[,order(col_data$RP, col_data$group)]  #ordering for supervised clustering on groups
  # make supervised heatmap
  pheatmap(e2, main = paste("Heatmap Supervised Clustering", gsub("_","",region)), scale = "row", show_rownames = F, show_colnames = F, cluster_cols=FALSE, cutree_rows = nc, annotation = annotation,
           annotation_row = pres.clust, annotation_names_row = F, annotation_colors = mycols, col = hmcol)
}

####################################################################
#                             CODE                                 #
####################################################################

# To make sure that the gene symbols of the Ensembl id's are known,
# a txt file is used of a gtf file (used during the alignment).
geneExp <- read.csv(file = "RawFiles/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)),
                    fill = T, na.strings=c("","NA"), stringsAsFactors = F)
# Only the gene id (ensemble id) and the gene name (genesymbol) was
# obtained from this file.
geneNames <- cbind(gsub("gene_id ", "", as.matrix(geneExp[,1])), gsub("gene_name ", "", as.matrix(geneExp[,5])),
                   trimws(apply(geneExp, 1, function(x) tail(na.omit(x), 1))))
gene_names <- as.data.frame(unique(geneNames[,1:2]))
colnames(gene_names) <- c("ensembl_gene_id", "external_gene_name")
# Make colors for heatmap use
mycols <- list(Group = c(C_HB_A = 'darkviolet', E4_HB_A = 'darkorchid', Ech_HB_A = 'deeppink', E1_HB_A = 'hotpink', C_HB_AG = 'navy',
                         E1_HB_AG = 'royalblue', E4_HB_AG = 'darkcyan', Ech_HB_AG = 'deepskyblue', C_SC_A = 'yellow', E1_SC_A = 'gold',
                         E4_SC_A = 'gold3', Ech_SC_A = 'orangered'),
               cluster = (c("cadetblue" , "hotpink", "navy", "orangered", "gold", "darkviolet")))
hmcol <- colorRampPalette(c("magenta", "dodgerblue2", "white", "tan3","firebrick3"))(n=99)

# Read input files and get nececary information
countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,c(4:6,13)]
colData$RP <- paste(colData$Region, colData$Population, sep = "_")

# start analyzer function for all file lists (all_results and each region)
count <- 0  #count to keep track of index in list
pdf(paste("Results/",prefix,"/",prefix,"_heatmaps.pdf",sep = ""))
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
  cat(paste("Start analysis for:", gsub("_","",region)))
  count = count + 1
  suppressWarnings(analyzer(files, region ))
  cat("\tdone!\n")
}
invisible(dev.off())

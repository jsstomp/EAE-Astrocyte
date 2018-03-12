####################################################################
# Author: Jafta Stomp
# Date: 15-02-2018)
# Description: 
#   This script finds differentially expressed genes using DESeq2
####################################################################
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least one argument must be supplied (prefix) (countfile).n", call.=FALSE)
} else if (length(args)==2) {
  prefix <- args[1]
  countfile <- args[2]
} else if (length(args)>2){
  stop("Too many arguments, please only supply (prefix) (countfile).n", call.=FALSE)
}

####################################################################
#                           IMPORTS                                #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))

# set working directory to parent parent working directory of this file
suppressMessages(setwd(gsub("DEA_DESeq2.R","",thisfile())))
setwd("../..")

# load count and col data
countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,4:6]
colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")

####################################################################
#                       DESeq2 analysis                            #
####################################################################

# Function to find differentially expressed genes between 2 given subgroups
do_DEA <- function(group1, group2) {
  columns <- colData[which(colData$Group == group1 | colData$Group == group2),]
  counts <- countData[,which(colnames(countData) %in% rownames(columns))]
  columns$Group <- as.factor(columns$Group)
  
  dds <- DESeqDataSetFromMatrix(counts,
                                columns,
                                design = ~ Group)
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.01, lfcThreshold = 2, pAdjustMethod = "BH")
  resLFC <- lfcShrink(dds, coef=2)
  de_genes <- res[which(res$log2FoldChange > 2 & res$padj < 0.05),]
  de_genes$ensembl_gene_id <- rownames(de_genes)
  file_name <- paste("DEA_Results/de_genes_", prefix, gsub("_", "", group1), "_vs_", gsub("_", "", group2), ".txt", sep = "")
  write.table(de_genes, file_name, quote=F, row.names=F, sep="\t")
}


#groups <- unique(colData$Group)
comparisons <- as.data.frame(t(data.frame(c("C_HB_A", "C_HB_AG"),c("C_HB_A", "C_SC_A"),c("C_SC_A", "C_HB_AG"),
                          c("E1_HB_A", "E1_HB_AG"),c("E1_HB_A", "E1_SC_A"),c("E1_SC_A", "E1_HB_AG"),
                          c("E4_HB_A", "E4_HB_AG"),c("E4_HB_A", "E4_SC_A"),c("E4_SC_A", "E4_HB_AG"),
                          c("Ech_HB_A", "Ech_HB_AG"),c("Ech_HB_A", "Ech_SC_A"),c("Ech_SC_A", "Ech_HB_AG"),
                          c("C_HB_A", "E1_HB_A"),c("E1_HB_A", "E4_HB_A"),c("E4_HB_A", "Ech_HB_A"),
                          c("C_HB_A", "E4_HB_A"),c("C_HB_A", "Ech_HB_A"),
                          c("E1_HB_A", "Ech_HB_A"),
                          c("C_HB_AG", "E1_HB_AG"),c("E1_HB_AG", "E4_HB_AG"),c("E4_HB_AG", "Ech_HB_AG"),
                          c("C_HB_AG", "E4_HB_AG"),c("C_HB_AG", "Ech_HB_AG"),
                          c("E1_HB_AG", "Ech_HB_AG"),
                          c("C_SC_A", "E1_SC_A"),c("E1_SC_A", "E4_SC_A"),c("E4_SC_A", "Ech_SC_A"),
                          c("C_SC_A", "E4_SC_A"),c("C_SC_A", "Ech_SC_A"),
                          c("E1_SC_A", "Ech_SC_A"), row.names =  c("group1","group2"))))
rownames(comparisons) <- 1:length(rownames(comparisons))
comparisons$group1 <- as.character(comparisons$group1)
comparisons$group2 <- as.character(comparisons$group2)

mapply(do_DEA, comparisons$group1, comparisons$group2)


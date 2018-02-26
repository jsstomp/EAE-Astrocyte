####################################################################
# Author: Jafta Stomp
# Date: 15-02-2018)
# Description: 
#   This script finds differentially expressed genes using DESeq2
####################################################################

setwd("~/Projects/EAE")
####################################################################
#                           IMPORTS                                #
####################################################################
library(DESeq2)
library(edgeR)
library(ggplot2)

countData <- read.table("Results/count_data_alt.txt", header = T, row.names = 1, check.names = F)
ph_data <- read.csv("RawFiles/target.csv", row.names = 4)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,4:6]
colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")

####################################################################
#                       DESeq2 analysis                            #
####################################################################
# Make DESeq object

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
  file_name <- paste("DEA_Results/de_genes_", gsub("_", "", group1), "_vs_", gsub("_", "", group2), ".txt", sep = "")
  write.table(de_genes, file_name, quote=F, col.names=NA, row.names=T, sep="\t")
}

do_DEA("C_HB_A", "C_HB_AG")
groups <- unique(colData$Group)
comparisons <- data.frame()
comparisons <- within(colData, {
  group1 = 
})
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


#!/usr/bin/env Rscript
####################################################################
# Author: Jafta Stomp
# Date: 28-06-2018
# Description: 
#   This script finds differentially expressed genes using DESeq2
####################################################################

####################################################################
#                           PARSER                                 #
####################################################################
suppressMessages(library(argparser))

parser <- arg_parser('R script for running DEA using DESeq2')
parser <- add_argument(parser, 'prefix', 
                       help='Please give a word or short description (no spaced) to be used in output files and directories')
parser <- add_argument(parser, 'countfile', help='absolute path to countfile')
parser <- add_argument(parser, 'FDR', type='numeric', help='wanted false discovery rate threshold')
parser <- add_argument(parser, 'logFC', type='numeric', help='wanted log fold change threshold')
parser <- add_argument(parser, '--analyze-results', flag = TRUE, 
                       help='optional argument if user wants to analyze results further')
parser <- add_argument(parser, '--venn', flag = TRUE, 
                       help='optional argument for if the user wants to create a venn diagram of interesting genes per region')
parser <- add_argument(parser, '--gene-ontology', flag = TRUE, 
                       help='optional argument for if the user wants to Gene Ontoloy Identification')
parser <- add_argument(parser, '--wgcna', flag = TRUE, 
                       help='optional argument for if the user wants to Weighted gene correlation network analysis (WGCNA)')

p <- parse_args(parser, argv=commandArgs(trailingOnly=TRUE))

prefix <- p$prefix[1]
countfile <- p$countfile[1]
FDR <- p$FDR[1]
logFC <- p$logFC[1]
analysis <- p$analyze_results[1]
venn <- p$venn[1]
GO <- p$gene_ontology[1]
WGCNA <- p$wgcna[1]


####################################################################
#                           IMPORTS                                #
####################################################################
cat("Loading required libraries.\n")
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))
suppressMessages(library(fpc))

# set working directory to parent parent working directory of this file
suppressMessages(setwd(gsub("DEA_DESeq2.R","",thisfile())))
setwd("../..")  #move out of git repository (where data is stored and saved to)

# make output directory
cat("Creating output directories.\n")
dir.create(paste("Results/", prefix, sep=""), showWarnings = FALSE)
dir.create(paste("Results/", prefix,"/DEA_Results/", sep=""), showWarnings = FALSE)

# load count and col data
cat("Reading input data.\n")
countData <- read.table(countfile, header = T, check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,4:6]
colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")

####################################################################
#                       DESeq2 analysis                            #
####################################################################
cat("Building functions.\n")
# Function to find differentially expressed genes between 2 given subgroups
do_DEA <- function(group1, group2) {
  columns <- colData[which(colData$Group == group1 | colData$Group == group2),]
  counts <- countData[,which(colnames(countData) %in% rownames(columns))]
  columns$Group <- as.factor(columns$Group)
  
  dds <- DESeqDataSetFromMatrix(counts,
                                columns,
                                design = ~ Group)
  dds <- DESeq(dds)
  res <- results(dds, alpha=FDR, lfcThreshold = logFC, pAdjustMethod = "BH")
  resLFC <- lfcShrink(dds, coef=2)
  de_genes <- res[which(abs(res$log2FoldChange) > logFC & res$padj < FDR),]
  de_genes$ensembl_gene_id <- rownames(de_genes)
  #distinguish up and down regulated genes
  up_genes <- res[which(res$log2FoldChange > logFC & res$padj < FDR),]
  up_genes$ensembl_gene_id <- rownames(up_genes)
  down_genes <- res[which(res$log2FoldChange < -logFC & res$padj < FDR),]
  down_genes$ensembl_gene_id <- rownames(down_genes)
  
  #file_name <- paste("Results/", prefix,"/DEA_Results", "/de_genes_", prefix, "_", gsub("_", "", group1), "_vs_", gsub("_", "", group2), ".txt", sep = "")
  file_name_up <- paste("Results/", prefix,"/DEA_Results", "/de_genes_", prefix, "_", gsub("_", "", group2), "_vs_", gsub("_", "", group1), ".txt", sep = "")
  file_name_down <- paste("Results/", prefix,"/DEA_Results", "/de_genes_", prefix, "_", gsub("_", "", group1), "_vs_", gsub("_", "", group2), ".txt", sep = "")
  #write.table(de_genes, file_name, quote=F, row.names=F, sep="\t")
  write.table(up_genes, file_name_up, quote=F, row.names=F, sep="\t")
  write.table(down_genes, file_name_down, quote=F, row.names=F, sep="\t")
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

cat(paste("Starting to do differential expression analysis on ",length(rownames(comparisons)),
          " different comparisons, please wait, this may take several minutes.\n",sep=""))
invisible(suppressMessages(mapply(do_DEA, comparisons$group1, comparisons$group2)))

if(analysis==TRUE){
  cat("Initiating DEA_analysis.R.\n")
  source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/DEA_analysis.R")
}
if(venn==TRUE){
  cat("Initiating venn_DEA.R.\n")
  source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/venn_DEA.R")
}
if(GO==TRUE){
  cat("Initiating GOI.R.\n")
  source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/GOI.R")
}
if(WGCNA==TRUE){
  cat("Initiating WGCNA.R.\n")
  source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/WGCNA.R")
}

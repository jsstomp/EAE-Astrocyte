####################################################################
# original author: M. Dubbelaar
# Date: 11-dec-2017
# Edited by: Jafta Stomp  (15-02-2018)
# Description: 
#   This script filters out lowly expressed genes using DAFS.R
#   and also creates PCA plots to look at if the different conditions
#   cluster together.
####################################################################
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (prefix).n", call.=FALSE)
} else if (length(args)==1) {
  prefix <- args[1]
}
#prefix <- "EAE2"

####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))

# set working directory to parent parent working directory of this file
setwd(gsub("prepare_expression.R","",thisfile()))
setwd("../..")

#source("Code/DAFS.R")

####################################################################
#                              Code                                #
####################################################################
# Load count data of all samples and all genes
countData <- read.table("RawFiles/mergedCounts.txt", header = T, row.names = 1, check.names = F)

#if(!file.exists(paste("Results/filtered_counts_", prefix, ".txt", sep=""))){
#  createDAFSfile(countData, prefix)
#}

#mat <- read.table(paste("Results/filtered_counts_", prefix, ".txt", sep=""), header = T, check.names = F)
#mat_mlcpm <- aveLogCPM(mat)

### Other filtering method, https://support.bioconductor.org/p/75914/ (Gordon Smith)
M <- median(colSums(countData)) * 1e-6 # median library size (millions)
idx <- rowSums( cpm(countData) >= 5/M ) >= ncol(countData)/2 
nr_features <- length(which(idx == TRUE))
countData <- countData[idx, ]
mlcpm <- aveLogCPM(countData)

# Information about the different samples is written into the target
# file, the step loads all of the necessary information.
target <- read.table("RawFiles//target.csv", header = T, sep=",")
ph_data <- read.table("RawFiles//180219_Sample\ info.csv", header = T, sep=",")
colnames(ph_data)[2] <- "Sample.ID.Seq"
target <- merge(target, ph_data[,c(2,4,6,10:12)], by="Sample.ID.Seq")

# load file containing a list of samples to be excluded based on quality control
exclude <- read.table("Results/excludes.txt")[,3]
# exclude these samples
target <- target[-which(target$Sample.ID.Marissa %in% exclude),]
# drop unused levels from factors in target
target <- droplevels(target)

countData <- countData[,which(colnames(countData) %in% target$Sample.ID.Marissa)]

source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/PCA_plotter.R")

# Write results
write.table(countData,"Results/count_data_NEW.txt",quote=F,col.names=NA,row.names=T, sep="\t")
write.table(target,"Results/col_data_NEW.txt",quote=F,col.names=T,row.names=T, sep="\t")

pdf(paste("Results/Density_", prefix, ".pdf", sep = ""))
ggplot() + aes(mlcpm) +
  geom_histogram(binwidth=0.2, colour = "tomato4", fill="tomato3") +
  labs(title="Filter Lowly Expressed Features", subtitle = "Densityplot of the average log CPM after filtering genes with less than 5 counts in half of the samples", x = "logCPM", y = "Frequency") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))
#ggplot() + aes(mat_mlcpm) +
 # geom_histogram(binwidth=0.2, colour = "tomato4", fill="tomato3") +
 #labs(title="Filter Lowly Expressed Features", subtitle = "Densityplot of the average log CPM after DAFS filtering", x = "logCPM", y = "Frequency") +
 #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))
dev.off()

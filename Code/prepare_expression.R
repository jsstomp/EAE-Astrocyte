####################################################################
# Author: Jafta Stomp
# Date: 27-02-2018
# Description: 
#   This script filters out lowly expressed genes using DAFS.R
#   and also creates PCA plots to look at if the different conditions
#   cluster together.
####################################################################
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least one argument must be supplied (prefix) (filtering method){DAFS|Smith}.n", call.=FALSE)
} else if (length(args)==2) {
  prefix <- args[1]
  method.filter <- toupper(args[2])
} else if (length(args)>2){
  stop("Too many arguments, please only supply (prefix) (filtering method){DAFS|Smith}.n", call.=FALSE)
}

if(method.filter=="DAFS"|method.filter=="SMITH"){
  cat("R script prepare_expression.R has been succesfully initiated.\n")
} else {
  stop("Filtering method should either be DAFS or Smith")
}

####################################################################
#                            IMPORTS                               #
####################################################################
cat("Importing necessary libraries...\n")
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))

# set working directory to parent parent working directory of this file
suppressMessages(setwd(gsub("prepare_expression.R","",thisfile())))
setwd("../..")

source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/DAFS.R")

####################################################################
#                              Code                                #
####################################################################
cat("Loading count data...\n")
# Load count data of all samples and all genes
countData <- read.table("RawFiles/mergedCounts.txt", header = T, row.names = 1, check.names = F)

# Use if statements to see which filtering method should be used
if(method.filter=="DAFS"){
  if(!file.exists(paste("Results/filtered_counts_", prefix, ".txt", sep=""))){
    createDAFSfile(countData, prefix)
  }
  test <- readline("What is the preferred filtering method, DAFS or Smith")
  countData <- read.table(paste("Results/filtered_counts_", prefix, ".txt", sep=""), header = T, check.names = F)
} else if(method.filter=="SMITH"){
  ### Other filtering method, https://support.bioconductor.org/p/75914/ (Gordon Smith)
  cat("Filtering count data...\n")
  M <- median(colSums(countData)) * 1e-6 # median library size (millions)
  idx <- rowSums( cpm(countData) >= 5/M ) >= ncol(countData)/2 
  nr_features <- length(which(idx == TRUE))
  countData <- countData[idx, ]
}
mlcpm <- aveLogCPM(countData)

# Information about the different samples is written into the target
# file, the step loads all of the necessary information.
cat("Loading phenotype data...\n")
target <- read.table("RawFiles//target.csv", header = T, sep=",")
ph_data <- read.table("RawFiles//180219_Sample\ info.csv", header = T, sep=",")
colnames(ph_data)[2] <- "Sample.ID.Seq"
target <- merge(target, ph_data[,c(2,4,6,10:12)], by="Sample.ID.Seq")

# load file containing a list of samples to be excluded based on quality control
cat("Excluding low quality samples...\n")
exclude <- read.table("Results/excludes.txt")[,3]
# exclude these samples
target <- target[-which(target$Sample.ID.Marissa %in% exclude),]
# drop unused levels from factors in target
target <- droplevels(target)
# remove remark from target, its no longer interesting and has empty values.
target <- target[-9]

countData <- countData[,which(colnames(countData) %in% target$Sample.ID.Marissa)]

PCA_check <- ""
while(PCA_check!="y"|PCA_check!="n"){
  con <- file("stdin")
  cat("Would you like to get some PCA plots? [y/n]: ")
  PCA_check <- readLines(con,1)
  close(con)
  if(PCA_check=="y"){
    oldw <- getOption("warn")
    options(warn = -1)
    cat("Initiating PCA_plotter.R from source...\n")
    suppressMessages(source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/PCA_plotter.R"))
    options(warn = oldw)
    break
  } else if(PCA_check=="n"){
    break
  }
}

target[,1] <- factor(toupper(as.character(target[,1])))
target$groupID <- paste(target$Condition, target$Region, target$Population, sep="_")
target <- target[match(colnames(countData), target$Sample.ID.Marissa),]

# Write results
cat("Writing results...\n")
write.table(countData,paste("Results/count_data_",prefix,".txt",sep=""),quote=F,col.names=NA,row.names=T, sep="\t")
write.table(target,"Results/col_data.txt",col.names=T,row.names=F, sep="\t")

# Density plot of filtered counts
cat("Create density plot of filtered counts.\n")
pdf(paste("Results/Density_", prefix, ".pdf", sep = ""))
ggplot() + aes(mlcpm) +
  geom_histogram(binwidth=0.2, colour = "tomato4", fill="tomato3") +
  labs(title="Filter Lowly Expressed Features", subtitle = "Densityplot of the average log CPM after filtering genes with less than 5 counts in half of the samples", x = "logCPM", y = "Frequency") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))
dev.off()

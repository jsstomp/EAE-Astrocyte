####################################################################
# Author: Jafta Stomp
# Date: 15-02-2018)
# Description: 
#   This script finds differentially expressed genes using DESeq2
####################################################################


####################################################################
#                           IMPORTS                                #
####################################################################
library(DESeq2)
library(ggplot2)

countData <- read.table("RawFiles/mergedCounts.txt", header = T, row.names = 1, check.names = F)
colData <- 

####################################################################
#                       DESeq2 analysis                            #
####################################################################
# Make DESeq object

dds <-DESeqDataSetFromMatrix()
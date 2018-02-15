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

####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))

# set working directory to parent parent working directory of this file
setwd(gsub("prepare_expression.R","",thisfile()))
setwd("..")

source("Code/DAFS.R")
####################################################################
#                    Load Necessary information                    #
####################################################################
# Load count data of all samples and all genes
countData <- read.table("RawFiles/mergedCounts.txt", header = T, row.names = 1, check.names = F)

createDAFSfile(countData, prefix)

mat <- read.table(paste("Results/filtered_counts_", prefix, ".txt", sep=""), header = T, check.names = F)

# Information about the different samples is written into the target
# file, the step loads all of the necessary information.
target <- read.table("RawFiles//target.csv", header = T, sep=",")

target[,1] <- factor(toupper(as.character(target[,1])))
target$groupID <- paste(target$Condition, target$Region, target$Population, sep="_")
target <- target[match(colnames(mat), target$Sample.ID.Marissa),]

####################################################################
#                          Preprocess Data                         #
####################################################################
target$groupID <- factor(target$groupID) 
dds <- DESeqDataSetFromMatrix(as.matrix(mat), target, ~groupID)
dds <- dds[ rowSums(counts(dds)) > 1, ]
####################################################################
#                             PCA plots                            #
####################################################################
# rlog transformation, may take some time, so please be patient.
#rld <- rlog(dds, blind=FALSE)
rld <- vst(dds)

data_ACSA <- plotPCA(rld, intgroup=c("groupID"), returnData=TRUE)
data_ACSA$Region <- target$Region
data_ACSA$Population <- target$Population
data_ACSA$Condition <- target$Condition
data_ACSA$SampleID <- target$Sample.ID.Seq
percentVar <- round(100 * attr(data_ACSA, "percentVar"))

cbPalette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffc948','#b15928')
colors <- scale_color_brewer(palette="Paired")
col <- colorRampPalette(brewer.pal(8, "Dark2"))(12)
hsv(3,2,1)

pdf(paste("Results/PCA_", prefix, ".pdf", sep = ""))
ggplot(data_ACSA, aes(PC1, PC2, color=Population)) +
  geom_point(size=5) +
  geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
  labs(title="PCA plot", subtitle="Population of astrocytes", x=paste0("PC1: ",percentVar[1],"% variance"), y=paste("PC2: ",percentVar[2],"% variance")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))
ggplot(data_ACSA, aes(PC1, PC2, color=Region)) +
  geom_point(size=5) +
  geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
  labs(title="PCA plot", subtitle="Region of astrocytes", x=paste0("PC1: ",percentVar[1],"% variance"), y=paste("PC2: ",percentVar[2],"% variance")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) 
ggplot(data_ACSA, aes(PC1, PC2, color=Condition)) +
  geom_point(size=5) +
  geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
  labs(title="PCA plot", subtitle="Condition of astrocytes", x=paste0("PC1: ",percentVar[1],"% variance"), y=paste("PC2: ",percentVar[2],"% variance")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) 
ggplot(data_ACSA, aes(PC1, PC2, color=groupID)) +
  geom_point(size=5) +
  geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
  labs(title="PCA plot", subtitle="Group of astrocytes", x=paste0("PC1: ",percentVar[1],"% variance"), y=paste("PC2: ",percentVar[2],"% variance")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  scale_color_manual(values=cbPalette)
dev.off()

####################################################################
# Author: Jafta Stomp
# Date: 12-03-2018
# Description: 
#   This script finds differentially expressed genes using DESeq2
#   And uses these results to make a venn diagram
#   only looks at control vs all other stages per region
####################################################################
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("At least one argument must be supplied (prefix) (countfile) (FDR) (logFC).n", call.=FALSE)
} else if (length(args)==4) {
  prefix <- args[1]
  countfile <- args[2]
  FDR <- args[3]
  logFC <- args[4]
} else if (length(args)>4){
  stop("Too many arguments, please only supply (prefix) (countfile) (FDR) (logFC).n", call.=FALSE)
}

####################################################################
#                           IMPORTS                                #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(kimisc))
suppressMessages(library(VennDiagram))
suppressMessages(library(gridExtra))

# set working directory to parent parent working directory of this file
suppressMessages(setwd(gsub("venn_DEA.R","",thisfile())))
setwd("../..")

# load count and col data
countData <- read.table(countfile, header = T, check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)
colnames(countData) <- rownames(ph_data)
colData <- ph_data[,4:6]
colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")

####################################################################
#                           FUNCTIONS                              #
####################################################################

venn_DEA <- function(group1, group2) {
  columns <- colData[which(colData$Group == group1 | colData$Group == group2),]
  counts <- countData[,which(colnames(countData) %in% rownames(columns))]
  columns$Group <- as.factor(columns$Group)
  
  dds <- DESeqDataSetFromMatrix(counts,
                                columns,
                                design = ~ Group)
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.01, lfcThreshold = 2, pAdjustMethod = "BH")
  return(res[c(2,6)])
}

vennDataBuilder <- function(db, region){
  # c-E1
  control = paste("C",region, sep="_")
  E1 = paste("E1", region, sep="_")
  E4 = paste("E4", region, sep="_")
  Ech = paste("Ech", region, sep="_")
  for(x in list(E1,E4,Ech)){
    pl = venn_DEA(control,x)
    db[paste("padj_",gsub(paste("_",region,sep=""),"",x),sep="")] <- pl$padj
    db[paste("log2FoldChange_",gsub(paste("_",region,sep=""),"",x),sep="")] <- pl$log2FoldChange
  }
  db = within(db, {
    E1.control = ifelse(log2FoldChange_E1 > logFC & padj_E1 < FDR, 1, 0)
    E4.control = ifelse(log2FoldChange_E4 > logFC & padj_E4 < FDR, 1, 0)
    Ech.control = ifelse(log2FoldChange_Ech > logFC & padj_Ech < FDR, 1, 0)
  })
  return(db)
}

venning <- function(db, title){
  a1 = nrow(subset(db, E1.control==1))
  a2 = nrow(subset(db, E4.control==1))
  a3 = nrow(subset(db, Ech.control==1))
  a12 = nrow(subset(db, E1.control==1 & E4.control==1))
  a13 = nrow(subset(db, E1.control==1 & Ech.control==1))
  a23 = nrow(subset(db, E4.control==1 & Ech.control==1))
  a123 = nrow(subset(db, E1.control==1 & E4.control==1 & Ech.control==1))
  
  g = draw.triple.venn(area1 = a1,area2 = a2, area3 = a3, n12 = a12, n23 = a23, n13 = a13, n123 = a123,
                   category = c("E1","E4","Ech"), fill = c("blue","red","green"), euler.d=T,scaled=T,ind=F)
  grid.arrange(gTree(children=g), top=title, bottom="VennDiagram of EAE stages vs control")
}

####################################################################
#                               CODE                               #
####################################################################

HBA_db <- data.frame(row.names=rownames(countData))
HBAG_db <- data.frame(row.names=rownames(countData))
SCA_db <- data.frame(row.names=rownames(countData))

suppressMessages(HBA_db <- vennDataBuilder(HBA_db,"HB_A"))
suppressMessages(HBAG_db <- vennDataBuilder(HBAG_db, "HB_AG"))
suppressMessages(SCA_db <- vennDataBuilder(SCA_db, "SC_A"))

pdf(paste("Results/", prefix, "_VennDiagrams.pdf",sep=""))
venning(HBA_db, "Region: Hindbrain acsa")
venning(HBAG_db, "Region: Hindbrain acsa-glast")
venning(SCA_db, "Region: Spinal cord acsa")
dev.off()

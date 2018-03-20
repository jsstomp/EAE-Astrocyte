####################################################################
# Author: Jafta Stomp
# Date: 12-03-2018
# Description: 
#   This script finds differentially expressed genes using DESeq2
#   And uses these results to make a venn diagram
#   only looks at control vs all other stages per region
####################################################################
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)<4) {
#   stop("At least one argument must be supplied (prefix) (countfile) (FDR) (logFC).n", call.=FALSE)
# } else if (length(args)==4) {
#   prefix <- args[1]
#   countfile <- args[2]
#   FDR <- as.numeric(args[3])
#   logFC <- as.numeric(args[4])
# } else if (length(args)>4){
#   stop("Too many arguments, please only supply (prefix) (countfile) (FDR) (logFC).n", call.=FALSE)
# }

####################################################################
#                           IMPORTS                                #
####################################################################
cat("Loading required libraries.\n")
suppressMessages(library(kimisc))
suppressMessages(library(VennDiagram))
suppressMessages(library(gridExtra))

# # set working directory to parent parent working directory of this file
# suppressMessages(setwd(gsub("venn_DEA.R","",thisfile())))
# setwd("../..")


####################################################################
#                           FUNCTIONS                              #
####################################################################

venn_de_reader <- function(group1, group2) {
  res <- read.table(paste("Results/", prefix, "/DEA_Results", "/de_genes_", prefix, "_", group1, "_vs_", group2, ".txt", sep=""), header = T, sep="\t")
  return(res[c(7)])
}

# vennDataBuilder is a function that calls on venn_DEA and compares the control of the given region with E1, E4 and Ech
# saving the FDR and logFC and adding them to the dataframe of the given region. Then checks within the dataframe which
# samples and genes have significantly de genes given the FDR and logFC thresholds given on the commandline.
# Returns the dataframe
vennDataBuilder <- function(db, region){
  # c-E1
  df <- data.frame()
  control = paste("C",region, sep="")
  E1 = paste("E1", region, sep="")
  E4 = paste("E4", region, sep="")
  Ech = paste("Ech", region, sep="")
  gene_count = 0
  for(x in list(E1,E4,Ech)){
    print(x)
    pl <- venn_de_reader(control,x)
    
    db[x] <- ifelse(rownames(db) %in% pl$ensembl_gene_id, 1, 0)
    
  }
  
  return(db)
}

# venning is a function that creates a venn diagram given a dataframe. The diagram consist of three circles (E1,E4,Ech)
venning <- function(db, region, title){
  E1.control = paste("E1", region, sep="")
  E4.control = paste("E4", region, sep="")
  Ech.control = paste("Ech", region, sep="")
  a1 = nrow(subset(db, get(E1.control)==1))
  a2 = nrow(subset(db, get(E4.control)==1))
  a3 = nrow(subset(db, get(Ech.control)==1))
  a12 = nrow(subset(db, get(E1.control)==1 & get(E4.control)==1))
  a13 = nrow(subset(db, get(E1.control)==1 & get(Ech.control)==1))
  a23 = nrow(subset(db, get(E4.control)==1 & get(Ech.control)==1))
  a123 = nrow(subset(db, get(E1.control)==1 & get(E4.control)==1 & get(Ech.control)==1))
  
  g = draw.triple.venn(area1 = a1,area2 = a2, area3 = a3, n12 = a12, n23 = a23, n13 = a13, n123 = a123,
                       category = c("E1","E4","Ech"), fill = c("blue","red","green"), euler.d=F, scaled=F, ind = F)
  grid.arrange(gTree(children=g), top=title, bottom="VennDiagram of EAE stages vs control")
}

####################################################################
#                               CODE                               #
####################################################################

# initiate the dataframes.
HBA_db <- data.frame(row.names=rownames(countData))
HBAG_db <- data.frame(row.names=rownames(countData))
SCA_db <- data.frame(row.names=rownames(countData))

# build dataframes using vennDataBuilder
cat("Running differential expression analysis for HB_A. please wait...\n")
suppressMessages(HBA_db <- vennDataBuilder(HBA_db,"HBA"))
cat("Running differential expression analysis for HB_AG. please wait...\n")
suppressMessages(HBAG_db <- vennDataBuilder(HBAG_db, "HBAG"))
cat("Running differential expression analysis for SC_A. please wait...\n")
suppressMessages(SCA_db <- vennDataBuilder(SCA_db, "SCA"))

# make vennDiagrams
cat("Creating Venn Diagrams and writing to Results.\n")
pdf(paste("Results/", prefix, "/", prefix, "_VennDiagrams.pdf",sep=""))
venning(HBA_db, "HBA", "Region: Hindbrain acsa")
venning(HBAG_db, "HBAG", "Region: Hindbrain acsa-glast")
venning(SCA_db,"SCA", "Region: Spinal cord acsa")
invisible(dev.off())
cat("Done!\n")

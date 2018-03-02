library(GSVA)
library(GSEABase)

#read analysis results
r1 <- read.table("DEA_Results/de_genes_E4SCA_vs_E4HBAG.txt")
set_name <- "E4SCA_vs_E4HBAG"

#create a geneset
writeLines(paste(set_name, paste(rownames(r1),collapse="\t"), sep="\t"), "r1.gmt")
geneSets <- getGmt("r1.gmt")

enrichment.scores <- gsva(data.matrix(countData), geneSets, method="gsva", mx.diff=TRUE, verbose=TRUE, parallel.sz=8)

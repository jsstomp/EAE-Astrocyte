####################################################################
# Author: Jafta Stomp
# Date: 26-03-2018
# Description: 
#   This script tries to identify Gene Ontology for a list of dea
#   results. It does this using a package topGO using the fisher 
#   exact test.
####################################################################
####################################################################
#                           IMPORTS                                #
####################################################################
suppressMessages(library(topGO))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(ggplot2))

####################################################################
#                           FUNCTIONS                              #
####################################################################
ens2eg <- function(ens_gene_list){
  cat("Converting ensembl gene IDs to entrez IDs...\t")
  eg.list <- as.list(org.Mm.egENSEMBL2EG)
  my.eg.list <- as.vector(unlist(eg.list[ens_gene_list]))
  cat("DONE\n")
  return(my.eg.list)
}

topGOing <- function(file_name, dea_name){
  t <- read.table(file_name, header=T,sep="\t",row.names=1)
  my.ensmusg <- ens2eg(t$ensembl_gene_id)
  geneList <- factor(as.integer (all_genes %in% my.ensmusg))
  names(geneList) <- all_genes
  cat("Building topGO object, this may take several minutes please wait...\t")
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "entrez")
  cat("DONE\n")
  sel.terms <- usedGO(GOdata)
  tt <- termStat(GOdata, sel.terms)
  tt <- tt[which(tt$Significant > 0),]
  tt <- tt[order(tt$Significant, decreasing=T),]
  df <- GOtest(GOdata)
  cat(paste("Writing results to ", paste("Results/", prefix,"/GO_Results/", sep=""), "...\t", sep=""))
  write.csv(df, file = paste("Results/", prefix,"/GO_Results/", dea_name, ".csv",sep=""))
  cat("DONE\n")
  return(df)
}

GOtest <- function(godat){
  #weight method
  cat("Running Fisher exact test, this may take several minutes please wait...\t")
  resultWeight <- runTest(godat, algorithm = "weight", statistic = "fisher")
  #classic method
  resultFisher <- runTest(godat, algorithm = "classic", statistic = "fisher")
  cat("DONE\n")
  pvalFis <- score(resultFisher)
  
  # #correlation between methods
  # pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
  # cor(pvalFis, pvalWeight)
  
  geneData(resultWeight)
  allRes <- GenTable(godat, classic = resultFisher, weight = resultWeight, 
                     orderBy = "weight", ranksOf = "classic", topNodes = 20)
  return(allRes)
}

plot_GO <- function(go_table, dea_name){
  # take -10log of p-value of fisher exact test for clearer visual
  go_table$min10log.weight <- -log10(as.numeric(go_table$weight))
  # Make sure table is ordered correctly
  go_table <- go_table[order(go_table$min10log.weight),]
  # factor to save order
  go_table$GO.ID <- factor(go_table$GO.ID, levels = go_table$GO.ID)
  go_table$Term <- factor(go_table$Term, levels = go_table$Term)
  # save levels to use for the labs
  reslabs1 <- levels(go_table$GO.ID)
  reslabs2 <- levels(go_table$Term)
  
  cat("Sending graphs to Results...")
  gp <- ggplot(go_table, aes(as.numeric(GO.ID),min10log.weight)) + #numeric GO.ID to use continues scale
    geom_bar(stat="identity", width=.5, fill="tomato3") +
    labs(title=dea_name,
         subtitle="Gene Ontology") + 
    coord_flip() +
    scale_y_continuous(name="-10log(p)") +
    scale_x_continuous(breaks=1:length(reslabs1), name = "GO",
                       labels=reslabs1,sec.axis=sec_axis(~.,
                                                         breaks=1:length(reslabs2),
                                                         labels=reslabs2)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
  png(paste("Results/", prefix,"/GO_Results/", dea_name, ".png",sep=""), width=800, height=800)
  print(gp)
  dev.off()
  cat("DONE\n")
}

# MAIN FUNCTION
main <- function(path){
  dea_name <- gsub(".txt","",gsub(paste("de_genes_", prefix, "_",sep=""),"",gsub(".*\\/\\/","",path)))
  f <- file(path, open="rb")
  # Check if nlines in file is > 1, if not, dont run functions
  nlines <- 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
    nlines <- nlines + sum(chunk == as.raw(10L))
  }
  invisible(close(f))
  
  if(nlines > 1){
    cat(paste("Starting Gene Enrichment Analysis on DEA results of: ", dea_name, ".\n", sep=""))
    df <- topGOing(path, dea_name)
    plot_GO(df, dea_name)
  }
  else{
    cat(paste(dea_name,"had no differentially expressed genes and will not be tested.\n"))
  }
}
####################################################################
#                               CODE                               #
####################################################################
#prefix <- "FDR001_logFC1"
# list files and read other files needed
file_list <- list.files(paste("Results/", prefix,"/DEA_Results/", sep=""), full.names=T)
countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
# create directory if it doesn't exist yet
dir.create(paste("Results/", prefix,"/GO_Results/", sep=""), showWarnings = FALSE)
# convert gene ensembl id's to entrez id's for all genes in countData
all_genes <- ens2eg(rownames(countData)) 

# Start main
suppressMessages(sapply(file_list,main))

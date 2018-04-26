####################################################################
# Author: Jafta Stomp
# Date: 25-04-2018
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
suppressMessages(library(annotate))

####################################################################
#                           FUNCTIONS                              #
####################################################################
ens2eg <- function(ens_gene_list){
  ## Function that converts ensemble gene ids to entrez gene ids
  ## Returns entrez gene ids
  
  cat("Converting ensembl gene IDs to entrez IDs...\t")
  eg.list <- as.list(org.Mm.egENSEMBL2EG)
  # map ensembl gene ids to a database of entrez ids
  my.eg.list <- as.vector(unlist(eg.list[ens_gene_list]))
  cat("DONE\n")
  return(my.eg.list)
}


topGOing <- function(file_name, dea_name){
  ## Function that builds a topGO object and retrieves gene symbols after testing
  ## Returns top 20 GOs for plotting
  
  dea_results <- read.table(file_name, header=T,sep="\t",row.names=1)
  # Convert ensmembl gene ids from dea results to entrez ids
  my.ezmusg <- ens2eg(dea_results$ensembl_gene_id)
  # Make a named list of genes and a 0 or 1 if gene occurs using all_genes
  geneList <- factor(as.integer (all_genes %in% my.ezmusg))
  names(geneList) <- all_genes
  cat("Building topGO object, this may take several minutes please wait...\t")
  # Build the topGO object, look only for Biological Process GOs in mouse
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "entrez")
  cat("DONE\n")
  # Initiate GOtest function to run fisher exact test
  df <- GOtest(GOdata, my.ezmusg)
  # Only keep results that have a weight of less than 0.01 and at least 3 significant genes mapped to it
  df <- df[which(df$weight < 0.01 & df$Significant >= 3),]
  # Convert entrez gene ids to HGNC gene names
  df$genes <- lapply(df$genes,function(genes){getSYMBOL(genes,data='org.Mm.eg')})
  # Collapse lists inside the columns and unlist the entire column to be able to write it to csv
  df$genes <- lapply(df$genes,function(genes){paste(genes,collapse=",")})
  df$genes <- unlist(df$genes)
  cat(paste("Writing results to ", paste("Results/", prefix,"/GO_Results/", sep=""), "...\t", sep=""))
  write.csv(df, file = paste("Results/", prefix,"/GO_Results/", dea_name, ".csv",sep=""))
  cat("DONE\n")
  # Only return the 20 most significant GOs to be plotted
  return(df[1:20,])
}

GOtest <- function(godat,significantGenes){
  ## Function that runs two types of fisher exact tests (weight and classic)
  ## Returns a dataframe with GO information
  
  # Weight method
  cat("Running Fisher exact test, this may take several minutes please wait...\t")
  resultWeight <- runTest(godat, algorithm = "weight", statistic = "fisher")
  # Classic method
  resultFisher <- runTest(godat, algorithm = "classic", statistic = "fisher")
  cat("DONE\n")
  
  # Build the data frame using GenTable function
  allRes <- GenTable(godat, classic = resultFisher, weight = resultWeight, 
                     orderBy = "weight", ranksOf = "classic", topNodes = 1000)
  # add only the genes that are in the DEA data (significantGenes)
  allRes$genes<-sapply(allRes$GO.ID, function(x) {
    genes<-genesInTerm(godat, x); genes[[1]][sapply(genes[[1]],function(id) id %in% significantGenes)]
  })
  
  return(allRes)
}

plot_GO <- function(go_table, dea_name){
  ## Function that plots the top 20 (or less) significant GOs in a horizontal barplot
  ## Prints the resulting plot to a result directory
  
  # take -10log of p-value of fisher exact test for clearer visual
  go_table$min10log.weight <- -log10(as.numeric(go_table$weight))
  # Make sure table is ordered correctly
  go_table <- go_table[order(go_table$min10log.weight),]
  # factor to save order
  go_table$GO.ID <- factor(go_table$GO.ID, levels = go_table$GO.ID)
  go_table$Term <- factor(go_table$Term, levels = go_table$Term)
  # save levels to use for the labs and breaks
  reslabs1 <- levels(go_table$GO.ID)
  reslabs2 <- levels(go_table$Term)
  cat("Sending graphs to Results...")
  gp <- ggplot(go_table, aes(as.numeric(GO.ID),min10log.weight)) + #numeric GO.ID to use continues scale
    geom_bar(stat="identity", width=.5, fill="tomato3") +
    labs(title=dea_name,
         subtitle="Gene Ontology") + 
    coord_flip() +
    scale_y_continuous(name="-10log(p)") +
    # Have GOid on left side of axis and GO term on right side using sex.acis
    scale_x_continuous(breaks=1:length(reslabs1), name = "GO",
                       labels=reslabs1,sec.axis=sec_axis(~.,
                                                         breaks=1:length(reslabs2),
                                                         labels=reslabs2)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),axis.text=element_text(size=12))
  png(paste("Results/", prefix,"/GO_Results/", dea_name, ".png",sep=""), width=800, height=800)
  print(gp)
  dev.off()
  cat("DONE\n")
}

# MAIN FUNCTION
main <- function(path){
  ## Main function of GOI.R, used to initiate other functions and set certain variables.
  ## Opens files if necesary.
  
  dea_name <- gsub(".txt","",gsub(paste("de_genes_", prefix, "_",sep=""),"",gsub(".*\\/\\/","",path)))
  # Open file to check the number of lines
  f <- file(path, open="rb")
  # Check if nlines in file is > 2, if not dont run functions (if file has only 1 gene, GO identification doesnt work)
  nlines <- 0L
  while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
    nlines <- nlines + sum(chunk == as.raw(10L))
  }
  invisible(close(f))
  
  if(nlines > 2){
    cat(paste("Starting Gene Enrichment Analysis on DEA results of: ", dea_name, ".\n", sep=""))
    df <- topGOing(path, dea_name)
    df <- df[!is.na(df),]
    if(!dim(df)[1] == 0){
      plot_GO(df, dea_name)
    }else cat("No significant results found, skipping the plotting step.\n")
  }
  else{
    cat(paste(dea_name,"had les than 2 differentially expressed genes and will not be tested.\n"))
  }
}
####################################################################
#                               CODE                               #
####################################################################
#prefix <- "FDR001_logFC1_all"
# list files and read other files needed
file_list <- list.files(paste("Results/", prefix,"/DEA_Results/", sep=""), full.names=T)
countData <- read.table("Results/count_data_smith.txt", header = T, check.names = F)
# create directory if it doesn't exist yet
dir.create(paste("Results/", prefix,"/GO_Results/", sep=""), showWarnings = FALSE)
# convert gene ensembl id's to entrez id's for all genes in countData
all_genes <- ens2eg(rownames(countData)) 

# Start main
suppressWarnings(suppressMessages(sapply(file_list,main)))

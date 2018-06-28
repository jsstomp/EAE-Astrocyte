####################################################################
# Author: Jafta Stomp
# Date: 28-02-2018
# Description: 
#   WGCNA script which runs WGCNA on expression data and returns 
#   plots of modules (clusters), as well as tables of top genes.
####################################################################
#                           IMPORTS                                #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(WGCNA))
suppressMessages(library(flashClust))

# import col and count data
smyth <- read.table("Results/count_data_smith.txt",header=T,check.names = F)
ph_data <- read.table("Results/col_data.txt", header = T, row.names = 4, check.names = F)


####################################################################
#                           FUNCTIONS                              #
####################################################################
main <- function(colData, countData){
  ddsL <- ddsBuilder(colData,countData)
  datExpr <- ddsL[[1]]
  colData <- ddsL[[2]]
  soft_power <- softThresholdPicker(datExpr)
  geneTree <- geneTreeBuilder(soft_power, datExpr)
  dynamicColors <- moduleMaker(geneTree)
  datTraits <- concensusAnalysis(colData,datExpr)
  datME <- eigenGeneAnalysis(datExpr,soft_power,dynamicColors)
  getModuleSignificance(datExpr,datTraits,dynamicColors)
  filteredGeneModules <- enrichmentAnalyzer(datExpr,dynamicColors)
  writeHubGenes(datExpr,datME,dynamicColors,filteredGeneModules)
}


ddsBuilder <- function(colData, countData){
  # Makes a dds object using col and countdata using DESeq2
  cat("Building DESeq object.\n")
  # filter and prepare col and count data for dds object
  ph_data <- colData[which(ph_data$Sample.ID.Marissa %in% colnames(countData)),]
  colnames(countData) <- rownames(ph_data)
  # rearrange columns to be ordered by region, population, condition
  order1 <- c(grep("C_HB_A_",colnames(countData)),grep("E1_HB_A_",colnames(countData)),grep("E4_HB_A_",colnames(countData)),grep("Ech_HB_A_",colnames(countData)),
              grep("C_HB_AG",colnames(countData)),grep("E1_HB_AG",colnames(countData)),grep("E4_HB_AG",colnames(countData)),grep("Ech_HB_AG",colnames(countData)),
              grep("C_S",colnames(countData)),grep("E1_S",colnames(countData)),grep("E4_S",colnames(countData)),grep("Ech_S",colnames(countData)))
  countData <- countData[,order1]
  
  colData <- ph_data[,4:6]
  colData <- colData[order1,]
  colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")
  countData <- countData[,which(colnames(countData) %in% rownames(colData))]
  
  ddsMat <- DESeqDataSetFromMatrix(countData,
                                   colData,
                                   design = ~ Group)
  dds<-DESeq(ddsMat)
  vsd <- getVarianceStabilizedData(dds)  
  
  vsd2 <- t(vsd)    
  datExpr<- vsd2
  
  return(list(datExpr,colData))
}


softThresholdPicker <- function(datExpr){
  # Automatically select soft power threshold using scale free topology fit
  cat("Selecting perfect soft threshold for WGCNA.\n")
  powers <- c(c(1:10), seq(from = 12, to=20, by=2));
  sft <- pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),
                           networkType = "signed")
  
  # plot results
  pdf(paste(output_directory,"power_plot.pdf",sep=""))
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
  
  # Red line corresponds to using an R^2 cut-off
  abline(h=0.8,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  # power chosen via sft$powerEstimate
  softPower <- sft$powerEstimate
  cat(paste("Chosen soft power threshold: ", softPower, ".\n",sep=""))
  return(softPower)
}


geneTreeBuilder <- function(softPower,datExpr){
  # builds similarity matrix and creates a genetree based on it
  cat("Building gene dendrogram.\n")
  # Makes an adjacency matrix based on genes (n*n), then creates a gene similarity matrix based on expression values
  # Returns a gene dendogram based on gene connectivity
  adj <- adjacency(datExpr,type = "signed",power=softPower)
  
  TOM <- TOMsimilarityFromExpr(datExpr = datExpr, networkType = "signed", power = softPower)
  
  colnames(TOM) =rownames(TOM) =colnames(datExpr)
  dissTOM <- 1-TOM
  
  geneTree <- flashClust(as.dist(dissTOM),method="average")
  
  pdf(paste(output_directory,"gene_tree_plot.pdf",sep=""))
  par(mfrow=c(1,1))
  plot(geneTree,xlab=",sub=",cex=0.3)
  dev.off()
  return(list(geneTree,dissTOM))
}


moduleMaker <- function(geneTree){
  # Creates modules based on the gene dendrogram
  cat("Searching for modules.\n")
  minModuleSize <- 200;
  gene_tree <- geneTree[[1]]
  dissTOM <- geneTree[[2]]
  # Module identification using dynamic tree cut
  dynamicColors <- labels2colors(cutreeDynamic(dendro = gene_tree,  method="tree", minClusterSize = minModuleSize, cutHeight = 0.95))
  staticColors <- as.character(cutreeStaticColor(dendro = gene_tree,cutHeight = 0.95, minSize = minModuleSize))
  hybridColors <- labels2colors(cutreeDynamic(dendro = gene_tree, distM = dissTOM, cutHeight = 0.95, deepSplit = 2, 
                                              pamRespectsDendro = FALSE, minClusterSize = minModuleSize))
  
  png(paste(output_directory,"gene_dendrogram_colors.png",sep=""), width = 1000, height = 1000, pointsize = 18, 
      type = "cairo")
  plotDendroAndColors(gene_tree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  sizeGrWindow(10,5)
  plotDendroAndColors(gene_tree, colors=data.frame(staticColors, dynamicColors, hybridColors), dendroLabels = FALSE,
                      marAll = c(1,8,3,1), main = "Gene dendrogram and module colors, TOM dissimilarity")
  dev.off()
  return(dynamicColors)
}


concensusAnalysis <- function(colData, datExpr){
  # performs consensus on the sample trait data by making a histogram of the sample expression levels
  cat("Applying consensus analysis.\n")
  sampleTree2 <- hclust(dist(datExpr), method = "average")
  colData$Group <- as.factor(colData$Group)
  
  colData.numeric <- apply(colData, 2, function(x) as.numeric(as.factor(x)))
  
  traitColors <- numbers2colors(colData.numeric,signed = F)
  
  pdf(paste(output_directory,"sample_dendro.pdf",sep=""))
  plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(colData),
                      main = "Sample dendrogram and trait heatmap", cex.dendroLabels = .6)
  dev.off()
  return(colData.numeric)
}


eigenGeneAnalysis <- function(datExpr, soft_power, dynamicColors){
  # Looks at the correlation between module eigengenes and initiates the creating of heatmaps.
  cat("Applying module eigengene analysis.\n")
  # Calculate the module eigengenes (1st PC) of modules in the dataset
  datME <- moduleEigengenes(datExpr,dynamicColors, softPower = soft_power)$eigengenes
  module_correlation <- signif(cor(datME, use="p"), 2)
  write.csv(module_correlation, paste(output_directory,"module_correlation.csv",sep=""))
  
  # We  define  a  dissimilarity  measure  between  the  module  eigengenes  that  keeps  track  of  the  sign  of  the  correlation
  # between the module eigengenes, and use it to cluster the eigengene:
  dissimME <- (1-t(cor(datME, method="p")))/2
  hclustdatME <- hclust(as.dist(dissimME), method="average")
  
  # Plot the eigengene dendrogram
  pdf(paste(output_directory,"module_eigengene_dendro.pdf",sep=""))
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based of the module eigengenes")
  dev.off()
  
  pdf(paste(output_directory,"module_histogram.pdf",sep=""))
  par(mfrow=c(1,1))
  plotMEpairs(datME)
  dev.off()
  
  png(paste(output_directory,"heatmap_modules_turq-gr-brwn.png",sep=""), width = 1000, height = 1000, pointsize = 28, 
      type = "cairo")
  par(mfrow=c(4,1), mar=c(1, 2, 7, 1), ps=16)
  multiHeatmapPlotter(c("blue","green","turquoise","red"),datExpr,dynamicColors)
  dev.off()
  png(paste(output_directory,"heatmap_modules_red_yel_blu.png",sep=""),, width = 1000, height = 1000, pointsize = 28, 
      type = "cairo")
  par(mfrow=c(4,1), mar=c(1, 2, 7, 1), ps=16)
  multiHeatmapPlotter(c("yellow","black","brown","green"),datExpr,dynamicColors)
  dev.off()
  
  singleHeatmapPlotter(dynamicColors,datExpr,datME)
  
  return(datME)
}


getModuleSignificance <- function(datExpr,datTraits.numeric,dynamicColors){
  # calculates the module significance of the modules based connectivity
  cat("Calculating average module significance.\n")
  # Plot module average gene significance per Trait (EAEscore,Region,Population and Group)
  pdf(paste(output_directory,"gene_significance.pdf",sep=""))
  par(mfrow = c(2,2),mar=c(5.1, 4.1, 4.1, 2.1))
  lapply(colnames(datTraits.numeric),function(coln){
    moduleSignificancePlotter(datExpr,datTraits.numeric,dynamicColors,coln)
  })
  dev.off()
}


moduleSignificancePlotter <- function(datExpr,datTraits.numeric,dynamicColors,trait){
  # plots module significance
  cat("Printing module significance plot.\n")
  GS <- as.numeric(cor(datExpr, factor(as.data.frame(datTraits.numeric)[[trait]]), use="p"))
  GeneSignificance <- abs(GS)
  plotModuleSignificance(GeneSignificance,dynamicColors, main= trait)
}


multiHeatmapPlotter <- function(four_modules,datExpr,dynamicColors){
  # plots heatmaps based on which modules were given as input
  cat("Printing heatmaps.\n")
  lapply(four_modules,function(which.module){
    plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
            clabels=colnames(t(scale(datExpr[,dynamicColors==which.module ]) )),rcols=which.module, cex.axis=1.5)
    title(which.module, line=6)
  })
}


singleHeatmapPlotter <- function(dynamicColors,datExpr,datME){
  # creates a combination of heatmap and barplot of a module
  cat("Printing more heatmaps.\n")
  lapply(names(table(dynamicColors)),function(module){
    ME=datME[, paste("ME",module, sep="")]
    png(paste(output_directory,"heat_bar_",module,".png",sep=""), width = 1000, height = 1000, pointsize = 18, 
        type = "cairo")
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 6, 2))
    plotMat(t(scale(datExpr[,dynamicColors==module ]) ),
            nrgcols=30,rlabels=F,rcols=module, clabels=colnames(t(scale(datExpr[,dynamicColors==module ]) )),
            cex.main=2)
    title(module, line=5)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME, col=module, main="", cex.main=2,
            ylab="eigengene expression",xlab="array sample")
    dev.off()
  })
}


enrichmentAnalyzer <- function(datExpr,dynamicColors){
  # find enrichment pathways for genes within a module
  cat("Analyzing user list enrichment.\n")
  genes <- colnames(datExpr)
  geneExp <- read.csv(file = "RawFiles/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)),
                      fill = T, na.strings=c("","NA"), stringsAsFactors = F)
  # Only the gene id (ensemble id) and the gene name (genesymbol) was
  # obtained from this file.
  geneNames <- cbind(gsub("gene_id ", "", as.matrix(geneExp[,1])), gsub("gene_name ", "", as.matrix(geneExp[,5])),
                     trimws(apply(geneExp, 1, function(x) tail(na.omit(x), 1))))
  gene_names <- as.data.frame(unique(geneNames[,1:2]))
  colnames(gene_names) <- c("ensembl_gene_id", "external_gene_name")
  
  geneModule <- data.frame(ensembl_gene_id = genes, module = dynamicColors)
  geneModuleM <- merge(gene_names,geneModule,by="ensembl_gene_id")
  geneModuleData <- geneModuleM[match(geneModule$ensembl_gene_id,geneModuleM$ensembl_gene_id),]
  filteredGeneModules <- na.omit(geneModuleData)
  uniqueFilteredGeneModules <- filteredGeneModules[ !duplicated(filteredGeneModules$external_gene_name), ] 
  # write a gene list ordered by module
  write.csv(uniqueFilteredGeneModules[order(uniqueFilteredGeneModules$module),],
            paste(output_directory,"gene_per_module.csv",sep=""))
  listGenes <- trimws(toupper(as.character(uniqueFilteredGeneModules[,2])))
  categories <- as.character(uniqueFilteredGeneModules[,3])
  
  enrichmentResults <- userListEnrichment(listGenes, categories,
                                          nameOut = paste(output_directory,"power8UserListEnrichment.csv",sep=""),
                                          useBrainLists = T,useBloodAtlases = T, useStemCellLists = T, useBrainRegionMarkers = T, 
                                          useImmunePathwayLists = T, usePalazzoloWang = T, omitCategories = "grey", outputGenes = T)
  return(filteredGeneModules)
}


writeHubGenes <- function(datExpr,datME,dynamicColors,filteredGeneModules){
  # write results to files
  cat("Saving hub genes.\n")
  datKME <- signedKME(datExpr, datME, outputColumnName="MM.")
  datKME <- datKME[-which(duplicated(filteredGeneModules$external_gene_name))]
  #top50Finder <- function(module){
  attach(datKME)
  df_list <- lapply(unique(dynamicColors),function(module){
    df <- assign(paste("datKME",module,sep="_"),data.frame(
      ensemble_gene_id = filteredGeneModules$ensembl_gene_id[order(get(paste("MM.",module,sep="")), decreasing = T)][1:50],
      gene_name = trimws(filteredGeneModules$external_gene_name[order(get(paste("MM.",module,sep="")), decreasing = T)][1:50]),
      gene_connectivity = round(get(paste("MM.",module,sep=""))[order(get(paste("MM.",module,sep="")), decreasing = T)][1:50], 3),
      module = module))
    return(df)
  })
  detach()
  
  lapply(df_list,function(df){
    write.csv(df,paste("~/Projects/EAE/Results/WGCNA_results/top50genes_",df$module[1], ".csv",sep=""),
              quote = F, row.names = F)
  })
  
  # get the actual hub genes
  hubs <- chooseTopHubInEachModule(datExpr = datExpr, dynamicColors,omitColors = "grey",power = 10, type = "signed")
  write.csv(as.data.frame(hubs), paste(output_directory,"hub_genes.csv",sep=""))
  
}
####################################################################
#                              CODE                                #
####################################################################\

output_directory <- "Results/WGCNA_results/"
main(ph_data,smyth)


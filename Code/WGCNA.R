library(DESeq2)
library(edgeR)
library(WGCNA)
library(flashClust)

de_genes_SCA <- read.csv("~/Projects/EAE/Results/FDR001_logFC1_all/de_genes_SCA")
de_genes_SCA <- de_genes_SCA$ensembl_gene_id
mergedCounts <- read.table("~/Projects/EAE/RawFiles/mergedCounts.txt",header = T,row.names=1)
smyth <- read.table("~/Projects/EAE/Results/count_data_smith.txt",header=T,check.names = F)

ph_data <- read.table("~/Projects/EAE/Results/col_data.txt", header = T, row.names = 4, check.names = F)

ph_data <- ph_data[which(ph_data$Sample.ID.Marissa %in% colnames(smyth)),]

colnames(smyth) <- rownames(ph_data)
colData <- ph_data[,4:6]
colData$Group <- paste(colData$Condition, colData$Region, colData$Population, sep="_")

# colData <- colData[which(colData$Region == "SC"),]
smyth <- smyth[,which(colnames(smyth) %in% rownames(colData))]

# countMatrix <- smyth[which(rownames(smyth) %in% de_genes_SCA),]

ddsMat <- DESeqDataSetFromMatrix(smyth,
                              colData,
                              design = ~ Group)
dds<-DESeq(ddsMat)
datExpr0<- assay(dds)

vsd <- getVarianceStabilizedData(dds)  

vsd2 <- t(vsd)    
datExpr<- vsd2

powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

#plot results
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# power of 5 doubled due to signed network
softPower <- 10

adj <- adjacency(datExpr,type = "signed",power=softPower)

TOM <- TOMsimilarityFromExpr(datExpr = datExpr, networkType = "signed", power = softPower)

colnames(TOM) =rownames(TOM) =rownames(smyth)
dissTOM <- 1-TOM

geneTree <- flashClust(as.dist(dissTOM),method="average")

par(mfrow=c(1,1))
plot(geneTree,xlab=",sub=",cex=0.3)

# Set the minimum module size
minModuleSize <- 200;

# Module identification using dynamic tree cut

dynamicMods <- cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize, cutHeight = 0.95);
#dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
staticColors <- as.character(cutreeStaticColor(dendro = geneTree,cutHeight = 0.95, minSize = minModuleSize))
hybridColors <- labels2colors(cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight = 0.95, deepSplit = 2, 
                                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize))
# table(dynamicColors)
table(staticColors)
table(hybridColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
sizeGrWindow(10,5)
plotDendroAndColors(geneTree, colors=data.frame(staticColors, dynamicColors, hybridColors), dendroLabels = FALSE,
                    marAll = c(1,8,3,1), main = "Gene dendrogram and module colors, TOM dissimilarity")

### CONSENSUS ANALYSIS ###
traitData <- read.csv("~/Projects/EAE/RawFiles/target.csv")
traitData$Group <- paste(traitData$Condition,traitData$Region,traitData$Population,sep="_")
dim(traitData)
names(traitData)

allTraits <- traitData[, -c(1,2,3,8,9)]

mySamples <- rownames(datExpr)
traitRows <- match(mySamples, allTraits$Unique.sample.ID)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

collectGarbage()

sampleTree2 <- hclust(dist(datExpr), method = "average")
datTraits$Group <- as.factor(datTraits$Group)

datTraits.numeric <- apply(datTraits, 2, function(x) as.numeric(as.factor(x)))

traitColors <- numbers2colors(datTraits.numeric,signed = F)

plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap", cex.dendroLabels = .6)
########### eigengene significance

# Calculate the module eigengenes (1st PC) of modules in the dataset
datME <- moduleEigengenes(datExpr,dynamicColors)$eigengenes
signif(cor(datME, use="p"), 2)

# We  define  a  dissimilarity  measure  between  the  module  eigengenes  that  keeps  track  of  the  sign  of  the  correlation
# between the module eigengenes, and use it to cluster the eigengene:
dissimME <- (1-t(cor(datME, method="p")))/2
hclustdatME <- hclust(as.dist(dissimME), method="average")

# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

# sizeGrWindow(8,9)
plotMEpairs(datME)#, y=ModuleEigengeneNetwork1$y)

# signif(cor(datME,ModuleEigengeneNetwork1[,-1]),2)

par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
######################################################################################
which.module="turquoise";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="green";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="brown";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

which.module="red";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="yellow";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="blue";
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
######################################################################################

for (module in names(table(dynamicColors))){
  ME=datME[, paste("ME",module, sep="")]
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[,dynamicColors==module ]) ),
          nrgcols=30,rlabels=F,rcols=module,
          main=module, cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
  }

signif(cor(y,datME, use="p"),2)
cor.test(y, datME$MEbrown)

p.values = corPvalueStudent(cor(y,datME, use="p"), nSamples = length(y))

# Calculating average module significance
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)

par(mfrow = c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
plotModuleSignificance(GeneSignificance,dynamicColors,2)

################USERLISTENRICHMENT####################
# Get dataframe of genes and their linked modules
genes <- dds@rowRanges@partitioning@NAMES

# To make sure that the gene symbols of the Ensembl id's are known,
# a txt file is used of a gtf file (used during the alignment).
geneExp <- read.csv(file = "RawFiles/Mus_musculus.GRCm38.85.txt", sep=";", header=F, col.names = paste0("V",seq_len(20)),
                    fill = T, na.strings=c("","NA"), stringsAsFactors = F)
# Only the gene id (ensemble id) and the gene name (genesymbol) was
# obtained from this file.
geneNames <- cbind(gsub("gene_id ", "", as.matrix(geneExp[,1])), gsub("gene_name ", "", as.matrix(geneExp[,5])),
                   trimws(apply(geneExp, 1, function(x) tail(na.omit(x), 1))))
gene_names <- as.data.frame(unique(geneNames[,1:2]))
colnames(gene_names) <- c("ensembl_gene_id", "external_gene_name")

geneModule <- data.frame(ensembl_gene_id = genes, module = dynamicColors)
geneModule <- merge(gene_names,geneModule,by="ensembl_gene_id")
geneModuleData <- geneModule[order(geneModule[,2]),1:3]
filteredGeneModules <- na.omit(geneModuleData)

uniqueFilteredGeneModules <- filteredGeneModules[ !duplicated(filteredGeneModules$external_gene_name), ]  #take first row within each id

listGenes <- trimws(toupper(as.character(uniqueFilteredGeneModules[,2])))

categories <- as.character(uniqueFilteredGeneModules[,3])

enrichmentResults <- userListEnrichment(listGenes, categories,nameOut = "testGenesUserListEnrichment.csv",useBrainLists = T,useBloodAtlases = T,
                   useStemCellLists = T, useBrainRegionMarkers = T, useImmunePathwayLists = T, usePalazzoloWang = T, omitCategories = "grey", outputGenes = T)

###################NETWORK VISUALISATION########################
# hubs <- chooseTopHubInEachModule(datExpr, dynamicColors, omitColors = "grey",
#                                  power=5, type = "signed")
# # The ensemble id will be converted into the gene name after that.
# dfHubs <- data.frame(ensembl_gene_id = as.character(hubs), module = as.character(names(hubs)))
# 
# hubs <- merge(dfHubs,gene_names,by="ensembl_gene_id")
# 
# # The idea of this part in the code it to make a file that can be
# # used in VisANT. Each of the genes that are avaiable in the 
# # M2 dataset are saved as probes.
# probes <- colnames(datExpr)
# 
# # Each module is used in this function.
# for (module in names(table(dynamicColors))){
#   # The right module is selected.
#   inModule <- dynamicColors==module
#   # The probes that are available for the module is saved as modProbes.
#   modProbes <- probes[inModule]
#   # The distances among the genes within a module are saved into modTOM
#   modTOM <- TOM[inModule, inModule]
#   dimnames(modTOM) <- list(modProbes, modProbes)
#   # The function exportNetworkToVisant makes sure that the available data 
#   # is saved into a txt file in a way that it can be used into VisANT.
#   # exportNetworkToVisANT(modTOM, file= paste("VisOutput/VisANTInput-", module, ".txt", sep=""),
#   # weighted = T, threshold = 0, probeToGene = data.frame(trimws(as.character(uniqueFilteredGeneModules[,1])), listGenes))
#   exportNetworkToCytoscape(modTOM, edgeFile= paste("VisOutput/CytoscapeInput-", module, "-edge.txt", sep=""),
#                            nodeFile= paste("VisOutput/CytoscapeInput-", module, "-node.txt", sep=""),
#                            weighted = T, threshold = 0.02, nodeAttr=module)
#   # softConnectifity calculated the connectivity of each node that are found 
#   # within one module. This data will be saved in to the dataframe IMConn
#   IMConn <- softConnectivity(datExpr[, modProbes])
#   # The top 30 genes within this dataset are saved.
#   top <- rank(-IMConn) <= 30
#   # And a file that is saved.
#   # print(head(modTOM[top,top]))
#   # exportNetworkToVisANT(modTOM[top, top], file= paste("VisOutput/VisANTInput-", module, "-top50.txt", sep=""),
#   #                       weighted = T, threshold = 0, probeToGene = data.frame(trimws(as.character(uniqueFilteredGeneModules[,1])), listGenes))
#   exportNetworkToCytoscape(modTOM[top, top], edgeFile= paste("VisOutput/CytoscapeInput-", module, "-edge-top50.txt", sep=""),
#                            nodeFile= paste("VisOutput/CytoscapeInput-", module, "-node-top50.txt", sep=""),
#                            weighted = T, threshold = 0.02, nodeAttr=module)
# }
# # Directory is set to save the environment since some steps like
# # calculation the TOM takes a lot of time.
# save.image()
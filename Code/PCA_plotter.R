####################################################################
# Author: Jafta Stomp
# Date: 27-02-2018
# Description: 
#   This script filters out lowly expressed genes using DAFS.R
#   and also creates PCA plots to look at if the different conditions
#   cluster together.
####################################################################

####################################################################
#                            IMPORTS                               #
####################################################################
cat("Importing necessary libraries...\n")
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

####################################################################
#                           FUNCTIONS                              #
####################################################################
cat("Building functions...\n")
cbPalette <- c('#1f78b4','#e31a1c','#33a02c','#696969','#b2df8a','#fdbf6f',
               '#ff7f00','#b15928','#a6cee3','#cab2d6','#fb9a99','#6a3d9a')

# Function to be called for building PCA plots of AEA astrocyte data.
# Uses PCdata and changes color depending on data input.
PCAplot <- function(data, PC1, PC2, color, subtitle, percent_var, colname){
  gg <- ggplot(data, aes(PC1, PC2, color=color)) +
    geom_point(size=5) +
    #geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
    labs(title="PCA plot", subtitle=subtitle, x=paste0("PC1: ",percent_var[1],"% variance"),
         y=paste("PC2: ",percent_var[2],"% variance")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "black")) + 
    scale_color_manual(values=cbPalette, name = colname)
  col <- colnames(data[which(colnames(data)==color)])
  if(colname=="RIN" | colname=="RNA_input"){
    gg <- gg + scale_colour_gradient(low = "red", high = "green")
  }
  else if(colname=="Score" | colname=="Input_rate") {
    gg <- gg + scale_colour_gradient(low = "green", high = "red")
  }
  gg
}

####################################################################
#                    Load Necessary information                    #
####################################################################
cat("Loading necessary information...\n")
# Make seperate target data per subtype for deeper PCA
target_HBA <- target[which(target$Region == "HB" & target$Population == "A"),]
target_HBAG <- target[which(target$Region == "HB" & target$Population == "AG"),]
target_SCA <- target[which(target$Region == "SC" & target$Population == "A"),]
# Also make seperate count matrices per subtype for deeper PCA
cd_HBA <- countData[,which(colnames(countData) %in% target_HBA$Sample.ID.Marissa)]
cd_HBAG <- countData[ ,which(colnames(countData) %in% target_HBAG$Sample.ID.Marissa)]
cd_SCA <- countData[,which(colnames(countData) %in% target_SCA$Sample.ID.Marissa)]

####################################################################
#                          Preprocess Data                         #
####################################################################
cat("Preprocessing data...\n")
target$groupID <- factor(target$groupID) 
dds <- DESeqDataSetFromMatrix(as.matrix(countData), target, ~groupID)
dds <- dds[ rowSums(counts(dds)) > 1, ]
# Make dds objects for the sub datasets as well.
dds_HBA <- DESeqDataSetFromMatrix(cd_HBA, target_HBA, ~groupID)
dds_HBA <- dds_HBA[ rowSums(counts(dds_HBA)) > 1, ]
dds_HBAG <- DESeqDataSetFromMatrix(cd_HBAG, target_HBAG, ~groupID)
dds_HBAG <- dds_HBAG[ rowSums(counts(dds_HBAG)) > 1, ]
dds_SCA <- DESeqDataSetFromMatrix(cd_SCA, target_SCA, ~groupID)
dds_SCA <- dds_SCA[ rowSums(counts(dds_SCA)) > 1, ]
####################################################################
#                             PCA plots                            #
####################################################################
# variance stabilizing transformation,  instead of rlog to save some time, and the quality stays high.
rld <- vst(dds)
vsd_HBA <- vst(dds_HBA)
vsd_HBAG <- vst(dds_HBAG)
vsd_SCA <- vst(dds_SCA)

# Create PCA data for the large dataset
data_ACSA <- plotPCA(rld, intgroup=c("groupID"), returnData=TRUE)
data_ACSA$Region <- target$Region
data_ACSA$Population <- target$Population
data_ACSA$Condition <- target$Condition
data_ACSA$SampleID <- target$Sample.ID.Seq
percentVar <- round(100 * attr(data_ACSA, "percentVar"))

# Initiate PCA data building for all sub datasets.
data_HBA <- plotPCA(vsd_HBA, intgroup=c("groupID"), returnData=TRUE)
data_HBAG <- plotPCA(vsd_HBAG, intgroup=c("groupID"), returnData=TRUE)
data_SCA <- plotPCA(vsd_SCA, intgroup=c("groupID"), returnData=TRUE)

count <- 0  #initiates count to be used in coming for loop
df_list <- list(HBA=data_HBA, HBAG=data_HBAG, SCA=data_SCA)
# Start building PCA data for all subsets using a for loop per subset in df_list
cat("Building PCA data structure...\n")
for(i in df_list) {
  count = count + 1  #add +1 to count to represent the current index of the df_list
  # Use assign to make percentVar and add collumns using within per dataset.
  assign(paste("percentVar_", names(df_list)[count], sep=""), round(100 * attr(get(paste("data_", names(df_list)[count], sep="")), "percentVar")))
  assign(paste("data_", names(df_list)[count], sep=""), within(i, {
    SampleID = get(paste("target_", names(df_list)[count], sep=""))$Sample.ID.Seq
    RIN = get(paste("target_", names(df_list)[count], sep=""))$RQN
    Score = get(paste("target_", names(df_list)[count], sep=""))$Score
    RNA_input = get(paste("target_", names(df_list)[count], sep=""))$Input.RNA..ng.
    Input_rate = get(paste("target_", names(df_list)[count], sep=""))$X..of.total.input
    Date = get(paste("target_", names(df_list)[count], sep=""))$Date
  }))
}

# library(genefilter)
# ntop <- 500
# rv <- rowVars(assay(rld))
# select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
# mat <- t( assay(rld)[select, ] )
# pca <- prcomp( mat )
# 
# 
# pca_rotated <- psych::principal(mat, rotate="varimax", nfactors=2, scores=TRUE)
# print(pca_rotated$scores[1:5,])  # Scores returned by principal()
# percentVarimax <- round(100 * attr(pca_rotated, "percentVar"))
# 
# vm <- varimax(mat)
# str(vm)
# 
# PCAplot(as.data.frame(vm$loadings[,1:2]),vm$loadings[1],vm$loadings[2], data_ACSA$groupID, "rotated", percentVar, "groupID")
# PCAplot(as.data.frame(pca_rotated$scores), pca_rotated$scores[,1], pca_rotated$scores[,2], data_ACSA$groupID, "rotated", percentVar, "groupID")
# PCAplot(data_ACSA, PC1, PC2, data_ACSA$groupID, "rotated", percentVar, "groupID")
# PCAplot()

pdf(paste("Results/PCA_", prefix, ".pdf", sep = ""))
# plot Population, Condition, Region and GroupID
cat("Creating plots Population, Condition, Region and GroupID...\n")
lapply(names(data_ACSA[c(6:8,4)]), function(x) {
  plots <- PCAplot(data_ACSA, PC1, PC2, data_ACSA[x][[1]], paste("Astrocytes by", x), percentVar, x)
  print(plots)
})

# plot per region_population
df_list <- list(data_HBA=data_HBA, data_HBAG=data_HBAG, data_SCA=data_SCA)
count <- 0
for(i in df_list){
  count = count + 1
  cat(paste("Creating plots for", gsub("data_","",names(df_list)[count]), "...\n", sep=""))
  plots <- lapply(colnames(i[c(4,10,9,8,7,6)]), function(x){
    PCAplot(i, PC1, PC2, i[x][[1]], paste(gsub("data_", "", names(df_list)[count]), "astrocytes by", x), unlist(lapply(paste("percentVar_",gsub("data_","", names(df_list)), sep=""), get)[count]), x)
    })
  print(plots)
}
invisible(dev.off())
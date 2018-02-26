####################################################################
# original author: M. Dubbelaar
# Date: 11-dec-2017
# Edited by: Jafta Stomp  (15-02-2018)
# Description: 
#   This script filters out lowly expressed genes using DAFS.R
#   and also creates PCA plots to look at if the different conditions
#   cluster together.
####################################################################

####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

####################################################################
#                           FUNCTIONS                              #
####################################################################
cbPalette <- c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#b15928',
               '#cab2d6','#1f78b4','#6a3d9a','#33a02c','#e31a1c','#696969')

PCAplot <- function(data, PC1, PC2, color, subtitle, percent_var, colname){
  gg <- ggplot(data, aes(PC1, PC2, color=color)) +
    geom_point(size=5) +
    geom_text(aes(label=SampleID),hjust=-0.5, vjust=-0.5, size=3, col="black") +
    labs(title="PCA plot", subtitle=subtitle, x=paste0("PC1: ",percent_var[1],"% variance"),
         y=paste("PC2: ",percent_var[2],"% variance")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(color = "black"))
  col <- colnames(data[which(colnames(data)==color)])
  if(colname=="groupID" | colname=="Date"){
    gg <- gg + scale_color_manual(values=cbPalette)
  }
  else if(colname=="RIN" | colname=="RNA_input"){
    gg <- gg + scale_colour_gradient(low = "red", high = "green")
  }
  else if(colname=="Score" | colname=="Input_rate") {
    gg <- gg + scale_colour_gradient(low = "green", high = "red")
  }
  try(gg)
}

####################################################################
#                    Load Necessary information                    #
####################################################################
target[,1] <- factor(toupper(as.character(target[,1])))
target$groupID <- paste(target$Condition, target$Region, target$Population, sep="_")
target <- target[match(colnames(countData), target$Sample.ID.Marissa),]
# Make seperate target data per subtype for deeper PCA
target_HBA <- target[which(target$Region == "HB" & target$Population == "A"),]
target_HBAG <- target[which(target$Region == "HB" & target$Population == "AG"),]
target_SCA <- target[which(target$Region == "SC" & target$Population == "A"),]
# Also make seperate count matrices per subtype for deeper PCA
cd_HBA <- countData[,which(colnames(countData) %in% target_HBA$Sample.ID.Marissa)]
cd_HBAG <- countData[,which(colnames(countData) %in% target_HBAG$Sample.ID.Marissa)]
cd_SCA <- countData[,which(colnames(countData) %in% target_SCA$Sample.ID.Marissa)]

####################################################################
#                          Preprocess Data                         #
####################################################################
target$groupID <- factor(target$groupID) 
dds <- DESeqDataSetFromMatrix(as.matrix(countData), target, ~groupID)
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds_HBA <- DESeqDataSetFromMatrix(cd_HBA, target_HBA, ~groupID)
dds_HBA <- dds_HBA[ rowSums(counts(dds_HBA)) > 1, ]
dds_HBAG <- DESeqDataSetFromMatrix(cd_HBAG, target_HBAG, ~groupID)
dds_HBAG <- dds_HBAG[ rowSums(counts(dds_HBAG)) > 1, ]
dds_SCA <- DESeqDataSetFromMatrix(cd_SCA, target_SCA, ~groupID)
dds_SCA <- dds_SCA[ rowSums(counts(dds_SCA)) > 1, ]
####################################################################
#                             PCA plots                            #
####################################################################
# rlog transformation, may take some time, so please be patient.
#rld <- rlog(dds, blind=FALSE)
rld <- vst(dds)
vsd_HBA <- vst(dds_HBA)
vsd_HBAG <- vst(dds_HBAG)
vsd_SCA <- vst(dds_SCA)

data_ACSA <- plotPCA(rld, intgroup=c("groupID"), returnData=TRUE)
data_ACSA$Region <- target$Region
data_ACSA$Population <- target$Population
data_ACSA$Condition <- target$Condition
data_ACSA$SampleID <- target$Sample.ID.Seq
percentVar <- round(100 * attr(data_ACSA, "percentVar"))

dataBuilder <- function(name){
  temp <- paste("data_", name, sep="")
  print(temp)
  assign(temp, "test") 
}

data_HBA <- plotPCA(vsd_HBA, intgroup=c("groupID"), returnData=TRUE)
data_HBAG <- plotPCA(vsd_HBAG, intgroup=c("groupID"), returnData=TRUE)
data_SCA <- plotPCA(vsd_SCA, intgroup=c("groupID"), returnData=TRUE)

count <- 0
df_list <- list(HBA=data_HBA, HBAG=data_HBAG, SCA=data_SCA)
for(i in df_list) {
  count = count + 1
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

# plot Population, Condition and Region
pdf(paste("Results/PCA_", prefix, ".pdf", sep = ""))
lapply(names(data_ACSA[c(6:8,4)]), function(x) {
  plots <- PCAplot(data_ACSA, PC1, PC2, data_ACSA[x][[1]], paste("Astrocytes by", x), percentVar, x)
  print(plots)
})
  
dev.off()
# plot per region_population

df_list <- list(data_HBA=data_HBA, data_HBAG=data_HBAG, data_SCA=data_SCA)

count <- 0
pdf(paste("Results/PCA_regPos_", prefix, ".pdf", sep = ""))
for(i in df_list){
  count = count + 1
  plots <- lapply(colnames(i[c(4,10,9,8,7,6)]), function(x){
    #print(class(i[x][[1]]))
    PCAplot(i, PC1, PC2, i[x][[1]], paste(gsub("data_", "", names(df_list)[count]), "astrocytes by", x), unlist(lapply(paste("percentVar_",gsub("data_","", names(df_list)), sep=""), get)[count]), x)
    })
  print(plots)
}
dev.off()


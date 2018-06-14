####################################################################
# Author: Jafta Stomp
# Date: 25-04-2018
# Description: 
#   This script groups GO's together towards their first common ancestor
#   and initiates circos_plotter.R to make a circos plot of genes and GOs
####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(GOSemSim))
suppressMessages(library(GOSim))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(source("Experimental-autoimmune-encephalomyelitis-Astrocyte-RNA-seq-analysis/Code/circos_plotter.R"))

####################################################################
#                            FUNCTIONS                             #
####################################################################
main <- function(file_list){
  ## Main function of this script, initializes other functions
  ## returns a dataframe object.
  
  cat("Find GOids.\n")
  GOids <- file_list_reader(file_list)[[1]]
  genes <- file_list_reader(file_list)[2]
  cat("build similarity matrix.\n")
  # getTermSim returns the pairwise similarities between GOterms in the form of a similarity matrix
  sim_mat <- getTermSim(GOids) # uses default method: 'Relevance'
  cat("build GO groups.\n")
  group_list <- lister(sim_mat)
  
  group_list <- bottom_adder(GOids,group_list)
  
  total <- go_gene_relinking(group_list,genes)

  groups <- go_grouper(total, file_list)
  return(total)
}


bottom_adder <- function(GOids,group_list){
  ## Adds GO terms to the data frame that did not have any significant common ancestry with other GO terms
  ## These all match to the highest possible GO: Biological Process (GO:0008150)
  ## returns the adjusted data frame
  
  for(GO in GOids){
    if(GO %ni% group_list$GO1){
      tdf <- data.frame(GO,GO,"GO:0008150", 0, 0)
      colnames(tdf) <- c("GO1","GO2","match","depth","top_level")
      group_list <- rbind(group_list,tdf)
    }
  }
  return(group_list)
}


file_list_reader <- function(file_list){
  ## Reads a list of files containing GO information
  ## Also checks if the files exist, since not all analysis have GO Results.
  ## Returns GO terms
  
  # Initialize a count in order to keep track of the previous file in the file list
  count <- 1
  for(f in file_list){
    if(!file.exists(f)){
      # remove file from the list if it doesn't exist
      file_list <- file_list[-count]
    }
    else{
      # Go to the next file
      count <- count + 1
    }
  }
  # Read the files and append them to form 1 data frame
  GOs <- lapply(file_list, read.csv, header=T)
  GOs <- do.call(rbind, GOs)
  GOs <- GOs[order(GOs$weight),]
  
  # Only GO IDs are interesting so only keep those and return them. also only keep the levels (not interested in duplicates)
  GOids <- levels(as.factor(GOs$GO.ID))
  # Genes are also interesting, here we can already discover the number of genes and the total number or fisher exact value
  genes <- GOs[,c(2,4,5,9,10)]
  return(list(GOids,genes))
}


lister <- function(sim_mat){
  ## Builds a data frame object containing GO terms with a similarity of larger than similarity
  ## Initializes mgetMinimumSubsumer2
  ## Returns the resulting dataframe to main
  
  df <- as.data.frame(matrix(nrow = 3))
  for(rn in rownames(sim_mat)){
    for(cn in colnames(sim_mat)){
      s <- sim_mat[rn,cn]
      # Similarity has to be above similarity (can be changed) and below 0.99 in order to avoid matching a GO to itself.
      if(s > similarity & s < 0.99){
        l <- c(rn,cn,s)
        # Add the matching GO's to the data frame
        df[rn] <- l 
      }
    }
  }
  # Initializes mgetMinimumSubsumer2 with the data frame
  robj <- mgetMinimumSubsumer2(df)
  return(robj)
}


mgetMinimumSubsumer2 <- function(df){
  ## Finds the first common ancestor of each pair in the data frame
  ## Then for each match found it wil go into deeper 'depth' untill a depth of 3 has been reached
  ## Initializes recursive_mgetMinumumSubsumer to find the following matches
  ## Initializes top_level_calculator to find the top level/depth of each match
  ## returns a data frame with paired GO's, matched GO, depth and top_level
  
  depth <- 1  # starting depth of the matches
  ndf <- as.data.frame(matrix(nrow=0,ncol=4))
  names(ndf) <- c("GO1","GO2","match","depth")
  # loop through the columns of the data frame the first time to get the first matches at a depth of 1
  for(col in colnames(df)[-1]){
    ms <- getMinimumSubsumer(df[1,col],df[2,col])
    tdf <- data.frame(df[1,col],df[2,col],ms,depth)
    names(tdf) <- c("GO1","GO2","match","depth")
    ndf <- rbind(ndf,tdf)
  }
  # initiate both old and new skiplist (old has 1 element for them to not be of the same length initially)
  skiplist <- c("empty")
  skiplist_new <- c()
  # while loop to go use getminimumsubsumer recursively (change 3 to something less hardcoded)
  while(depth < 3){
    depth <- depth + 1
    # Break out of loop if no new matches are made before a depth of 3 has been reached
    if(length(skiplist_new)==length(skiplist)){
      break()
    }
    # Remember the last skiplist
    skiplist <- skiplist_new
    # Initiate recursive_mgetMinimumSubsumer
    ndf <- recursive_mgetMinimumSubsumer(ndf, depth, skiplist)
    # Remember the next skiplist
    skiplist_new <- ndf$GO1
    # Set max_depth to current depth
    max_depth <- depth
  }
  
  # Re-initiate depth at 1
  depth <- 1
  # Create a data_frame for saving the top levels
  level_list <- as.data.frame(matrix(nrow=0,ncol=2))
  # Until depth reaches max_depth calculate the top level for each match
  while(depth <= max_depth){
    # Initiate top_level_calculator for current depth
    top_levels <- top_level_calculator(ndf, depth)
    level_list <- rbind(level_list,top_levels)
    # Go to the next depth
    depth <- depth + 1
  }
  # Add the top levels to the original data frame
  ndf$top_level <- level_list$top_level
  
  # Some matches are still missing the correct top level higher up in the data frame
  done_list <- c()  # Checks if the matches have already been fixed
  for(match in ndf$match){
    if(match %in% ndf[which(ndf$depth==3),]$match & match %ni% done_list){
      if(ndf[which(ndf$match==match),]$top_level == 2){
        # Matches that have a top_level of 2 but a top level of 3 as well in lower matches get set to 3.
        ndf[which(ndf$match==match & ndf$top_level==2),]$top_level <-3
      }
      # Update done list
      done_list <- c(done_list,match)
    }
  }
  return(ndf)
}


top_level_calculator <- function(df, depth){
  ## Calculates if a match is used to find a new match as well.
  ## If not the current depth is the top_level
  ## Otherwise the top_level is raised by 1
  ## Returns a dataframe of top levels
  
  # Create a dataframe that can hold the top levels
  top_levels <- as.data.frame(matrix(nrow=0,ncol=2))
  to <- df[which(df$depth==depth),]
  to_plus <- df[which(df$depth==depth+1),]
  for(match in to$match){
    if(match %ni% to_plus$GO1 & match %ni% to_plus$GO2){
      # If match is not part of a new pair current depth is top level
      top_level <- depth
      tdf <- data.frame(match,top_level)
      top_levels <- rbind(top_levels,tdf)
    }
    else{
      # Else top level is depth + 1
      top_level <- depth + 1
      tdf <- data.frame(match, top_level)
      top_levels <- rbind(top_levels,tdf)
    }
  }
  return(top_levels)
}


recursive_mgetMinimumSubsumer <- function(df, depth, skiplist){
  ## Checks the previously made matches to see if they can also match eachother to form greater groups.
  ## Returns a dataframe 1 depth lower
  
  # Create a dataframe similar to input dataframe
  ndf <- as.data.frame(matrix(nrow=1,ncol=4))
  names(ndf) <- c("GO1","GO2","match","depth")
  for(go1 in df[,"match"]){
    for(go2 in df[,"match"]){
      s <- getTermSim(c(go1,go2))[1,2]
      # Check if GO's have not been used before and are not the same and have a similarity of more than similarity
      if(go1 != go2 & s > similarity & go1 %ni% skiplist & go2 %ni% skiplist){
        # Get the first common ancestor of the GO's
        ms <- getMinimumSubsumer(go1,go2)
        # Save the information of the GO's in a temporary data frame
        tdf <- data.frame(go1,go2,ms,depth)
        names(tdf) <- c("GO1","GO2","match","depth")
        if(!(tdf$GO1 %in% ndf$GO1 | tdf$GO1 %in% ndf$GO2) & !(tdf$GO2 %in% ndf$GO1 | tdf$GO2 %in% ndf$GO2)){
          # If the match has not yet been added to the data frame do so.
          ndf <- rbind(ndf,tdf)
        }
      }
    }
    # Add used GO to the skiplist
    skiplist <- c(skiplist, go1)
  }
  # merge the 2 dataframes and return
  df <- rbind(df,ndf[-1,])
  return(df)
}


go_gene_relinking <- function(gos,genes){
  colnames(gos)[1] <- "GO.ID"
  df <- merge(gos,genes,by="GO.ID")
  # print(dim(gos))
  # print(str(genes))
  return(df)
}


gene_deduplication <- function(df){
  nl <- list()
  member_list <- c()
  for(i in 1:length(rownames(df))){
    GOid <- as.character(df[i,]$GO.ID)
    genes <- unlist(strsplit(as.character(df[i,]$genes), ","))
    if(GOid %in% member_list){
      nl[[GOid]] <- unique(c(nl[[GOid]],genes))
    }
    else{
      nl[[GOid]] <- c(genes)
      member_list <- c(member_list,GOid)
    }
  }
  # Replace significant n genes in df by total significant unique genes
  for(x in df$GO.ID){
    df[which(df$GO.ID == x),]$Significant <- length(nl[[x]])
  }
  df$proportion.genes <- apply(df,1,function(x){
    round(as.numeric(x[7])/as.numeric(x[6]),digits = 3)
  })
}


go_grouper <- function(g, file_list){
  if(length(g[which(g$depth==3)]$match)<2){
    max_depth <- 2
  }else max_depth <- 3
  d3Gos <- c(levels(as.factor(as.character(g[which(g$depth==max_depth),]$match))),"GO:0008150")

  # data.frame object with 11 columns, first column is the different go_groups and the other 10 are the top 10 GO's matched to the group
  # NAs are initiated for when a group has less than 10 GO's matched
  pr <- get.genes(file_list)
  
  # determine the direction of the current set of files.
  direction <- NULL
  # if the analysis starts with control, it contains down regulated genes/ GOs
  if(startsWith(pr[[2]][1], "C")){
    direction <- "down"
  } else direction <- "up"
  
  go_groups <- data.frame(GOgroup = d3Gos, GO1 = NA, GO2 = NA, GO3 = NA, GO4 = NA, GO5 = NA,
                          GO6 = NA, GO7 = NA, GO8 = NA, GO9 = NA, GO10 = NA, ngenes = 0, genes = NA)
  # for every group try to find the GOs linked to it
  for(go in go_groups$GOgroup){
    all_gos <- g[which(g$match==go),c(1,2,8)]
    for(i in 1:length(rownames(all_gos))){
      all_gos <- rbind(all_gos,g[which(g$match==as.character(all_gos[i,1])),c(1,2,8)])
      all_gos <- rbind(all_gos,g[which(g$match==as.character(all_gos[i,2])),c(1,2,8)])
    }
    # Get number of genes for each group
    genes <- g[which(g$match==go),9]
    genes2 <- NULL
    for(genel in genes){
      genes2 <- c(genes2,strsplit(genel,","))
    }
    # Set column for number of genes
    go_groups[which(go_groups$GOgroup == go),12] <- length(unique(unlist(genes2)))
    # Set column for genes
    go_groups[[which(go_groups$GOgroup == go),13]] <- list(unique(unlist(genes2)))
    # for the circos plot, we want a seperate list of or df of only the top 10 most significant GOs per group
    top_list <- NULL
    # Therefor order all_gos on fisher exact weight value and select the top 10
    for(i in unique(all_gos$GO.ID)){
      s <- all_gos[which(all_gos$GO.ID == i),]
      r1 <- s[order(s$weight),][1,c(1,3)]
      r2 <- s[order(s$weight),][1,c(2,3)]
      colnames(r2) <- colnames(r1)
      top_list <- rbind(top_list, r1, r2)
    }
    o <- unique(top_list[order(top_list$weight),]$GO.ID)
    for(i in 1:length(o)){
      if(i > 10) break()
      else go_groups[which(go_groups$GOgroup == go),i+1] <- as.character(o[i])
    }
    
  }
  # the 3 EAE scores that are compared in the circosplot to control
  tests <- c("c1","c4","cch")
  # Get the number of genes per GO group for each test and total (list of data frames)
  go_group.list <- lapply(1:length(rownames(go_groups)),function(i){
    c1 <- counting(go_groups[i,13], pr[[1]][1])
    c4 <- counting(go_groups[i,13],pr[[1]][2])
    cch <- counting(go_groups[i,13],pr[[1]][3])
    df2 <- data.frame(go_groups = go_groups[i,1], c1 = c1, c4 = c4, cch = cch, total = c1+c4+cch)
  })
  # Bind the dataframes together by row
  rdf <- do.call(rbind,go_group.list)
  # Get total for each test as wel (colSums)
  cs <- colSums(rdf[,2:4])
  rdf[length(rownames(rdf))+1,] <- NA
  for(i in names(cs)){
    rdf[length(rownames(rdf)),i] <- cs[i]
  }
  # Make a proportional data frame to be used for making the circos plot proportions.
  tdf = NULL
  for(type in tests) {
    for(s in d3Gos) {
      l = rdf[[1]] == s
      n = sum(l)
      n = 1
      dd = data.frame(type = rep(type, n),
                      species = rep(s, n),
                      vaule1 = rdf[which(rdf$go_groups==s), type]/rdf[which(rdf$go_groups==s),"total"],
                      value2 = rdf[which(rdf$go_groups==s),type]/rdf[length(rownames(rdf)),type])
      tdf = rbind(tdf, dd)
    }
  }
  # Initialize circos_plotter.R and write to pdf.
  pdf(file = paste("Results/",prefix,"/",direction,"circos.pdf",sep=""))
  circler(g, d3Gos, tdf, go_groups, tests)
  dev.off()
  
  return(go_groups)
}


counting <- function(genes, l){
  # Function for counting the number of genes matching between GO group and test
  genes <- unlist(genes)
  l <- unlist(l)
  count <- 0
  for(i in genes){
    if(i %in% l){
      count <- count + 1
    }
  }
  return(count)
}


get.genes <- function(file_list){
  #open dea files to find all degs
  #convert deg from ensemble to gene names
  dea.prefixs <- unlist(lapply(file_list,function(i)lapply(strsplit(i,"/"), function(u)gsub(".csv","",u[4]))))
  dea.files <- paste("Results/",prefix,"/DEA_Results/de_genes_",prefix,"_",dea.prefixs,".txt",sep="")
  dea <- lapply(dea.files, read.delim, sep="\t", header=T)
  test <- lapply(dea,ens2symbol)
  return(list(test,dea.prefixs))
}

ens2symbol <- function(dea){
  # Function that converts ensembl gene ids to entrez gene ids and then to hgnc gene symbols
  # returns the gene names
  ens.genes <- dea$ensembl_gene_id
  enz.genes <- ens2eg(ens.genes)
  names.genes <- lapply(enz.genes,function(g){getSYMBOL(g,data='org.Mm.eg')})
  names.genes <- lapply(names.genes,function(g){paste(g,collapse=",")})
  names.genes <- unlist(names.genes)
  names.genes <- unique(names.genes)
  
  return(names.genes)
}

goIdToTerm <- function(x, names = TRUE, keepNA = TRUE) {
  stopifnot(requireNamespace("GO.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ans <- rep(NA_character_, length(x))
  names(ans) <- x
  ids <- AnnotationDbi::GOID(GO.db::GOTERM)
  i <- match(x, ids)
  k <- which(!is.na(i))
  res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
  ans[k] <- res
  if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
  if (!names) names(ans) <- NULL
  return(ans)
}


tableFormatter <- function(groups){
  tier2 <- as.character(unique(groups[which(groups$depth==2),]$match))
  cols_df <- c(tier2, goIdToTerm(tier2))
  
  df <- as.data.frame(matrix(nrow=100,ncol=length(tier2)*2))
  colnames(df) <- cols_df
  for(go in tier2){
    tier1 <- as.character(unlist(groups[which(groups$match==go),2:3]))
    tier0 <- NULL
    for(got1 in tier1){
      t0 <- as.character(unlist(groups[which(groups$match==got1),2:3]))
      tier0 <- c(tier0, t0)
    }
    tier0 <- unique(tier0)
    for(i in 1:length(tier0)){
      df[i,go] <- tier0[i]
      df[i,goIdToTerm(go)] <- goIdToTerm(tier0[i])
    }
    
  }
  df <- df[,unlist(lapply(1:length(tier2),function(i)seq(i,length(tier2)*2,length(tier2))))]
  return(df)
}
####################################################################
#                              CODE                                #
####################################################################
'%ni%' <- Negate('%in%')  # Easy to use reverse of the %in% function (not in)

# prefix <- "FDR005_logFC1_all"
go_results <- paste("Results/",prefix,"/GO_Results/",sep="")
group_results <- gsub("GO_Results/","",go_results)
similarity <- 0.6  # Current threshold for GO term similarity

# A list containing lists of files that contain GOs that need to be grouped together
filelist <- list(
  c(paste(go_results,c("CSCA_vs_E1SCA.csv","CSCA_vs_E4SCA.csv","CSCA_vs_EchSCA.csv"), sep="")),
  c(paste(go_results,c("E1SCA_vs_CSCA.csv","E4SCA_vs_CSCA.csv","EchSCA_vs_CSCA.csv"), sep=""))
)

# Mouse database
mmGO <- godata('org.Mm.eg.db', ont="BP")

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

#  

# Run main function for each list of files
groups <- lapply(filelist, function(fl){
  main(fl)
})

# Write files
for(i in 1:2){
  if(i==1){
    direct <- "down"
  }else if(i==2){
    direct <- "up"
  } 
  write.csv(groups[i], paste(group_results, direct, "_groups.csv",sep=""))
}

# reformat tables to human readable
up_groups <- read.csv(paste(group_results,"up_groups.csv",sep=""))
down_groups <- read.csv(paste(group_results,"down_groups.csv",sep=""))

up_df <- tableFormatter(up_groups)
down_df <- tableFormatter(down_groups)

# write again
write.csv(up_df, paste(group_results,"up_groups_table.csv",sep=""))
write.csv(down_df, paste(group_results,"down_groups_table.csv",sep=""))



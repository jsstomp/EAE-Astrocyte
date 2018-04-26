####################################################################
# Author: Jafta Stomp
# Date: 25-04-2018
# Description: 
#   This script groups GO's together towards their first common ancestor
####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(GOSemSim))
suppressMessages(library(GOSim))

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
  genes <- GOs[,c(2,4,5,10)]
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
  # while loop to go use getminimumsubsumer recursively (change 10 to something less hardcoded)
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
  print(dim(gos))
  print(str(genes))
  return(df)
}


gene_deduplication <- function(df){
  nl <- list()
  member_list <- c()
  for(i in 1:length(rownames(df))){
    GOid <- as.character(df[i,]$GO.ID)
    # print(as.character(GOid))
    genes <- unlist(strsplit(as.character(df[i,]$genes), ","))
    # print(genes)
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

####################################################################
#                              CODE                                #
####################################################################
'%ni%' <- Negate('%in%')  # Easy to use reverse of the %in% function (not in)

go_results <- "Results/FDR001_logFC1_all/GO_Results/"
similarity <- 0.6  # Current threshold for GO term similarity

# A list containing lists of files that contain GOs that need to be grouped together
filelist <- list(
  c(paste(go_results,c("CSCA_vs_E1SCA.csv","CSCA_vs_E4SCA.csv","CSCA_vs_EchSCA.csv"), sep="")),
  c(paste(go_results,c("E1SCA_vs_CSCA.csv","E4SCA_vs_CSCA.csv","EchSCA_vs_CSCA.csv"), sep=""))
)

# Mouse database
mmGO <- godata('org.Mm.eg.db', ont="BP")

# g <- main(filelist[[1]])

# Run main function for each list of files
groups <- lapply(filelist, function(fl){
  # print(fl)
  main(fl)
  # print(levels(main(fl)$GO_group_id))
})
levels(as.factor(groups[[1]][["GO_group_id"]])) #of down regulated genes
levels(as.factor(groups[[2]][["GO_group_id"]])) #of up regulated genes
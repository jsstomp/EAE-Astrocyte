####################################################################
# Author: Jafta Stomp
# Date: 29-03-2018
# Description: 
#   This script groups GO's together towards their parents
####################################################################
#                            IMPORTS                               #
####################################################################
suppressMessages(library(GOSemSim))
suppressMessages(library(GOSim))

####################################################################
#                            FUNCTIONS                             #
####################################################################
file_list_reader <- function(file_list){
  count <- 1
  for(f in file_list){
    if(!file.exists(f)){
      file_list <- file_list[-count]
    }
    else{
      count <- count + 1
    }
  }
  GOs <- lapply(file_list, read.csv, header=T)
  GOs <- do.call(rbind, GOs)
  GOs <- GOs[order(GOs$weight),]
  
  GOids <- levels(as.factor(GOs$GO.ID))
  return(GOids)
}
# 
# build_sim_matrix <- function(GOids){
#   df <- data.frame(matrix(nrow = 0, ncol= 3))
#   colnames(df) <- c("GO.1", "GO.2", "sem.sim")
#   # create initial groups
#   for(go1 in GOids){
#     for(go2 in GOids){
#       s <- goSim(go1,go2, mmGO)
#       de <- data.frame(go1,go2,s)
#       colnames(de) <- c("GO.1", "GO.2", "sem.sim")
#       df <- rbind(df,de)
#     }
#   }
#   sim_mat <- with(df, tapply(sem.sim, list(GO.1,GO.2), mean))
#   print(sim_mat)
#   break()
#   return(sim_mat)
# }

is.contained <- function(vec1,vec2){
  x=vector(length = length(vec1))
  for (i in 1:length(vec1)) {
    x[i] = vec1[i] %in% vec2
    if(length(which(vec1[i] %in% vec2)) == 0) vec2 else
      vec2=vec2[-match(vec1[i], vec2)]
  }
  y=all(x==T)
  return(y)
}

lister <- function(sim_mat){
  # print(sim_mat)
  # list_list <- list()
  # for(coln in colnames(sim_mat)){
  #   l <- list(coln)
  #   row_count <- 1
  #   l_count <- 2
  #   for(s in sim_mat[,coln]){
  #     if(s > 0.6 & s < 0.99){
  #       l[l_count] <- rownames(sim_mat)[row_count]
  #       l_count <- l_count + 1
  #     }
  #     row_count <- row_count + 1
  #   }
  #   list_list[[coln]] <- l
  # }
  df <- as.data.frame(matrix(nrow = 3))
  for(rn in rownames(sim_mat)){
    for(cn in colnames(sim_mat)){
      s <- sim_mat[rn,cn]
      if(s > 0.6 & s < 0.99){
        l <- c(rn,cn,s)
        df[rn] <- l
      }
    }
  }
  robj <- mgetMinimumSubsumer2(df)
  # print(robj)
}

mgetMinimumSubsumer2 <- function(df){
  depth <- 1
  ndf <- as.data.frame(matrix(nrow=0,ncol=4))
  names(ndf) <- c("GO1","GO2","match","depth")
  for(col in colnames(df)[-1]){
    ms <- getMinimumSubsumer(df[1,col],df[2,col])
    tdf <- data.frame(df[1,col],df[2,col],ms,depth)
    names(tdf) <- c("GO1","GO2","match","depth")
    ndf <- rbind(ndf,tdf)
  }
  # depth <- 2
  # skiplist <- c()
  # len <- length(rownames(ndf))
  # count <- 0
  # for(go1 in ndf[-1,"match"]){
  #   for(go2 in ndf[-1,"match"]){
  #     s <- getTermSim(c(go1,go2))[1,2]
  #     if(go1 != go2 & s > 0.6 & go2 %ni% skiplist){
  #       ms <- getMinimumSubsumer(go1,go2)
  #       tdf <- data.frame(go1,go2,ms,depth)
  #       names(tdf) <- c("GO1","GO2","match","depth")
  #       ndf <- rbind(ndf,tdf)
  #     }
  #   }
  #   skiplist <- c(skiplist, go1)
  #   count <- count + 1
  #   if(count == len){
  #     depth <- depth + 1
  #     len <- length(rownames(ndf)) -len
  #     count <- 0
  #   }
  # }
  
  # initiate both old and new skiplist (old has 1 element for them to not be of the same length initially)
  skiplist <- c("empty")
  skiplist_new <- c()
  # while loop to go use getminimumsubsumer recursively (change 10 to something less hardcoded)
  # ndf <- ndf[-1,]
  while(depth < 10){
    depth <- depth + 1
    # Break out of loop if no new matches are made.
    if(length(skiplist_new)==length(skiplist)){
      break()
    }
    skiplist <- skiplist_new
    ndf <- recursive_mgetMinimumSubsumer(ndf, depth, skiplist)
    skiplist_new <- ndf$GO1
    max_depth <- depth - 1
  }

  depth <- 1
  level_list <- as.data.frame(matrix(nrow=0,ncol=2))
  while(depth <= max_depth){
    top_levels <- TBN(ndf, depth)
    level_list <- rbind(level_list,top_levels)
    depth <- depth + 1
  }
  # print(level_list)
  ndf <- merge(ndf,level_list,by=match)
  #ndf$top_level <- level_list
  print(ndf)
  return(ndf)
}

TBN <- function(df, depth){
  top_levels <- as.data.frame(matrix(nrow=0,ncol=2))
  memory_list <- c()
  print(depth)
  for(match in df[which(df$depth == depth),]$match){
    # print(df[which(df$match == match),]$depth)
    # print(match)
    # if(df[which(df$match == match),]$depth == depth){
    if(match %ni% df$GO1 | match %ni% df$GO2){
      top_level <- depth
      tdf <- data.frame(match,top_level)
      top_levels <- rbind(top_levels,tdf)
    }
    else{
      memory_list <- c(memory_list,match)
    }
    # }
  }
  return(top_levels)
}

recursive_mgetMinimumSubsumer <- function(df, depth, skiplist){
  # o_skiplist <- c()
  ndf <- as.data.frame(matrix(nrow=1,ncol=4))
  names(ndf) <- c("GO1","GO2","match","depth")
  for(go1 in df[,"match"]){
    skip2 <- "go2"
    for(go2 in df[,"match"]){
      s <- getTermSim(c(go1,go2))[1,2]
      if(go1 != go2 & s > 0.6 & go1 %ni% skiplist & go2 %ni% skiplist){
        ms <- getMinimumSubsumer(go1,go2)
        tdf <- data.frame(go1,go2,ms,depth)
        names(tdf) <- c("GO1","GO2","match","depth")
        ndf <- rbind(ndf,tdf)
      }
      skip2 <- go2
    }
    skiplist <- c(skiplist, go1)
    
  }
  df <- rbind(df,ndf[-1,])
  return(df)
}

'%ni%' <- Negate('%in%')

mgetMinimumSubsumer <- function(sim_mat){

  if(length(go_group) == 1){
    name <- go_group[1]
  }
  else if(length(go_group) == 2){
    name <- getMinimumSubsumer(go_group[1],go_group[2])
  }
  else{
    l <- list()
    for(go1 in go_group){
      for(go2 in go_group){
        if(!go1==go2){
          ms<-getMinimumSubsumer(go1,go2)
          l[[ms]] <- ms
        }
      }
    }
    name <- mgetMinimumSubsumer(unlist(l))
  }
  # test <- as.data.frame(go_group)
  # print(test)
  
  # test <- lapply(go_group, function(g){
  #   c(g,name)
  # })
  # print(test)
  # cat("all GO in group\n")
  # print(name)
  # cat("GO-group\n")
  # print(go_group)
  return(name)
}

write_results <- function(df){
  
  
}
main <- function(file_list){
  # print(file_list)
  cat("Find GOids.\n")
  GOids <- file_list_reader(file_list)
  
  cat("build similarity matrix.\n")
  # sim_mat <- build_sim_matrix(GOids)
  sim_mat <- getTermSim(GOids[1:100]) # Remove 1:100
  
  cat("build GO groups.\n")
  group_list <- lister(sim_mat)
  # print(length(group_list))
  # print(str(group_list))
  # df <- data.frame(GOids)
  # print(group_list)
  # print(group_list)
  break()
  go_groups<-as.data.frame(t(as.data.frame(lapply(group_list, function(l){
    # print(str(group_list))
    if(length(l)>2){
      m <- mgetMinimumSubsumer(l)
      # print(m)
      return(m)
    }
    else{
      # print("other")
      # return as highest GO (Biological Process)
      return("GO:0008150")
    }
    
  }))))
  colnames(go_groups) <- "GO_group_id"
  #levels(go_groups$GO_group_id)
  fname <- gsub(go_results,"",gsub("HBAG.csv","",gsub("HBAG_vs","",file_list[1])))
  # print(fname)
  fname <- paste(go_results,paste("grouped_GO_",fname,".csv",sep=""),sep="")
  # print(fname)
  write.csv(go_groups, file = fname)
  # print(levels(go_groups$GO_group_id))
  print(levels(go_groups$GO_group_id))
  # return(go_groups)
}
####################################################################
#                              CODE                                #
####################################################################
go_results <- "Results/FDR001_logFC1_all/GO_Results/"
# filelist <- list(
#   c(paste(go_results,c("CHBAG_vs_E1HBAG.csv","CHBA_vs_E1HBA.csv","CSCA_vs_E1SCA.csv"), sep="")),
#   c(paste(go_results,c("CHBAG_vs_E4HBAG.csv","CHBA_vs_E4HBA.csv","CSCA_vs_E4SCA.csv"), sep="")),
#   c(paste(go_results,c("CHBAG_vs_EchHBAG.csv","CHBA_vs_EchHBA.csv","CSCA_vs_EchSCA.csv"), sep="")),
#   c(paste(go_results,c("E1HBAG_vs_E4HBAG.csv","E1HBA_vs_E4HBA.csv","E1SCA_vs_E4SCA.csv"), sep="")),
#   c(paste(go_results,c("E1HBAG_vs_EchHBAG.csv","E1HBA_vs_EchHBA.csv","E1SCA_vs_EchSCA.csv"), sep="")),
#   c(paste(go_results,c("E4HBAG_vs_EchHBAG.csv","E4HBA_vs_EchHBA.csv","E4SCA_vs_EchSCA.csv"), sep=""))
# )
filelist <- list(
  c(paste(go_results,c("CSCA_vs_E1SCA.csv","CSCA_vs_E4HBA.csv","CSCA_vs_EchSCA.csv"), sep="")),
  c(paste(go_results,c("E1SCA_vs_CSCA.csv","E4SCA_vs_CSCA.csv","EchSCA_vs_CSCA.csv"), sep=""))
)
# filelist<-c(paste(go_results,c("CSCA_vs_E1SCA.csv","CSCA_vs_E4SCA.csv","CSCA_vs_EchSCA.csv",
                     # "E1SCA_vs_CSCA.csv","E4SCA_vs_CSCA.csv","EchSCA_vs_CSCA.csv"), sep=""))
# groups <- main(filelist)

# runs <- c("C_E1", "C_E4", "C_Ech", "E1-E4", "E1-Ech", "E4-Ech")
# filelist <- lapply(runs, function(r){
#   f1 <- gsub("_.*","", r)
#   f2 <- gsub(".*_","", r)
#   c(paste(go_results,c(paste(f1,"HBAG_vs_",f2,"HBAG.csv",sep=""),
#                      paste(f1,"HBA_vs_",f2,"HBA.csv",sep=""),
#                      paste(f1,"SCA_vs_",f2,"SCA.csv",sep="")),sep=""))
# })

mmGO <- godata('org.Mm.eg.db', ont="BP")



groups <- lapply(filelist, function(fl){
  # print(fl)
  main(fl)
  # print(levels(main(fl)$GO_group_id))
})
levels(as.factor(groups[[1]][["GO_group_id"]])) #of down regulated genes
levels(as.factor(groups[[2]][["GO_group_id"]])) #of up regulated genes

#################TESTING BELOW############################

# for each group, check which GO's (per item in file list) are in it.

list1 <- list()
for(g in groups){
  for(go in levels(g$GO_group_id)){
    
    g<-getTermSim(c(GOids[1], go))
    # g<- goSim(GOids[1], go, mmGO)
    if(g[1,2] > 0){
      list1 <- list(list1, list(g[1,2], go))
      print(go)
      print(g[1,2])
    }
    
  }
}
for(i in unlist(list1)){
  print(i)
}
print(max(unlist(list1)))


goSim("GO:0010745", "GO:0010875", mmGO)
getTermSim(c("GO:0010745", "GO:0010875"))
mat<-getTermSim(GOids, method="relevance")

apply(mat,c(1,2), function(s){
  if(s > 0.6 & s < 0.99){
    print(s)
  }
})

mapply(function(x,y){
  if(mat[x,y] > 0.6 & mat[x,y] < 0.99){
    l <- c(x,y,mat[x,y])
    print(l)
  }
}, rownames(mat), colnames(mat))

for(rown in rownames(mat)){
  for(coln in colnames(mat)){
    if(mat[rown,coln] > 0.6 & mat[rown,coln] < 0.99){
      l <- c(rown,coln,mat[rown,coln])
      print(l)
    }
  }
}
# GO:0050896

library(igraph)
# for(i in filelist){
#   GOids <- file_list_reader(i)
#   g <- getGOGraph(GOids)
#   g2 <- igraph.from.graphNEL(g)
#   plot(g2)
# }

# Nice for making visualization
GOids <- file_list_reader(filelist[[1]])
g <- getGOGraph(c("GO:0010745", "GO:0010875"))
g2 <- igraph.from.graphNEL(g)
GOids[c(1,3)]
plot(g2, vertex.label=V(g2)$name)
getMinimumSubsumer("GO:0010745", "GO:0035855")

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
    # print(file_list[[count]])
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
  
  GOids <- levels(GOs$GO.ID)
  return(GOids)
}

build_sim_matrix <- function(GOids){
  df <- data.frame(matrix(nrow = 0, ncol= 3))
  colnames(df) <- c("GO.1", "GO.2", "sem.sim")
  # create initial groups
  for(go1 in GOids){
    for(go2 in GOids){
      s <- goSim(go1,go2, mmGO)
      de <- data.frame(go1,go2,s)
      colnames(de) <- c("GO.1", "GO.2", "sem.sim")
      df <- rbind(df,de)
    }
  }
  sim_mat <- with(df, tapply(sem.sim, list(GO.1,GO.2), mean))
  return(sim_mat)
}

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
  list_list <- list()
  for(coln in colnames(sim_mat)){
    l <- list(coln)
    row_count <- 1
    l_count <- 2
    for(s in sim_mat[,coln]){
      if(s > 0.3 & s < 1){
        l[l_count] <- rownames(sim_mat)[row_count]
        l_count <- l_count + 1
      }
      row_count <- row_count + 1
    }
    list_list[[coln]] <- l
  }
  
  list_list <- list_list[order(sapply(list_list, length), decreasing=T)]
  
  list_list2 <- list()
  count <- 1
  skip_list <- c()
  for(a in list_list){
    a <- unlist(a)
    nlist <- list()
    if(!a[1] %in% skip_list){
      for(b in list_list){
        if(!b[1] %in% skip_list){
          b <- unlist(b)
          c <- is.contained(b,a)
          if(c==TRUE){
            nlist[[a[1]]] <- a 
            skip_list <- c(skip_list, b[1])
            count <- count + 1
          }
        }
      }
    }
    list_list2[[a[1]]] <- nlist
  }
  list_list2 <- list_list2[order(sapply(list_list2, length), decreasing=T)]
  
  list_list2 <- lapply(list_list2, unlist)
}

mgetMinimumSubsumer <- function(go_group){
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
  GOids <- file_list_reader(file_list)
  sim_mat <- build_sim_matrix(GOids)
  group_list <- lister(sim_mat)
  # print(length(group_list))
  # print(str(group_list))
  # df <- data.frame(GOids)
  # print(group_list)
  go_groups<-as.data.frame(t(as.data.frame(lapply(group_list, function(l){
    # print(str(group_list))
    if(length(l)>2){
      m <- mgetMinimumSubsumer(l)
      # print(m)
      return(m)
    }
    else{
      # print("other")
      return("other")
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
  return(go_groups)
}
main(filelist[[1]])
####################################################################
#                              CODE                                #
####################################################################
go_results <- "Results/FDR001_logFC1/GO_Results/"
filelist <- list(
  c(paste(go_results,c("CHBAG_vs_E1HBAG.csv","CHBA_vs_E1HBA.csv","CSCA_vs_E1SCA.csv"), sep="")),
  c(paste(go_results,c("CHBAG_vs_E4HBAG.csv","CHBA_vs_E4HBA.csv","CSCA_vs_E4SCA.csv"), sep="")),
  c(paste(go_results,c("CHBAG_vs_EchHBAG.csv","CHBA_vs_EchHBA.csv","CSCA_vs_EchSCA.csv"), sep="")),
  c(paste(go_results,c("E1HBAG_vs_E4HBAG.csv","E1HBA_vs_E4HBA.csv","E1SCA_vs_E4SCA.csv"), sep="")),
  c(paste(go_results,c("E1HBAG_vs_EchHBAG.csv","E1HBA_vs_EchHBA.csv","E1SCA_vs_EchSCA.csv"), sep="")),
  c(paste(go_results,c("E4HBAG_vs_EchHBAG.csv","E4HBA_vs_EchHBA.csv","E4SCA_vs_EchSCA.csv"), sep=""))
)
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

#################TESTING BELOW############################

# for each group, check which GO's (per item in file list) are in it.

for(g in groups){
  for(go in levels(g$GO_group_id)){
    print(go)
    g<-getTermSim(c(GOids[1], go))
    # g<- goSim(GOids[1], go, mmGO)
    print(g[1,2])
  }
}


getTermSim(c(GOids[1],"GO:0071704"))
GO:0050896

# library(igraph)
# for(i in filelist){
#   GOids <- file_list_reader(i)
#   g <- getGOGraph(GOids)
#   g2 <- igraph.from.graphNEL(g)
#   plot(g2)
# }

# Nice for making visualization
GOids <- file_list_reader(filelist[[1]])
g <- getGOGraph(c(GOids[1], "GO:0071704"))
g2 <- igraph.from.graphNEL(g)
GOids[c(1,3)]
plot(g2, vertex.label=V(g2)$name)
getMinimumSubsumer(GOids[1],"GO:0050896")

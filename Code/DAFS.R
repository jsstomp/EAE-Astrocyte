##############################################################################################
# qcut incorporated from: supplementory data of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4098771/
# createDAFSfile function by: Marissa Dubbelaar
# Edited by: Jafta Stomp (15-02-2018)
# Description:
#   This script is partialy incorporated from the official DAFS paper
#   qcut uses the one sample ks-test and multivariate adaptive regression splines to
#   estimate the optimal cutoff for filtering lowly expressed genes.
#   These cutoffs are used per gene and 10% of the samples per gene need to pass to be kept.
##############################################################################################
suppressMessages(library(mclust))  ###used to estimate the parameters of reference distribution 
suppressMessages(library(earth))   ###used to apply MARS algorithm
###############################################
#                  Functions                  #
###############################################
## Function that perform DAFS transformation
qcut <- function(data,name) {
  #set vector for cutoff values
  cutv <- rep(0,0)
  
  for(i in 1:ncol(data)) {
    xx <- data[,i]
    xx <- xx[-which(xx==0)]
    
    #take log2 of the data
    log2xx <- log2(xx)
    dlog2 <- data.frame(LogC=log2xx)
    
    #vector to store Kolmogorov Smirnov stats
    vv <- rep(0,0)
    
    #select start point
    start <- length(log2xx[log2xx==min(log2xx)])/length(log2xx)
    
    #set sequence
    s <- seq(round(start,2),0.5,by=0.005)
    
    #loop through cuts of the data to determine targeted K-S statistic
    for(q in s) {
      #select data greater than a quantile and run Mclust on that data to determine theoretical distro
      d <- log2xx[which(log2xx>quantile(log2xx,q,na.rm=T))]
      out <- Mclust(d,G=1)
      ks <- ks.test(d,"pnorm",out$parameter$mean,out$parameter$variance$sigmasq)
      vv <- c(vv,ks$statistic)
    }
    
    #determine first left-most local minima
    out <- earth(s,vv,thresh=0.005)
    
    #save suggested cut
    cutv <- c(cutv,min(out$cuts[out$cuts>0]))
  }
  
  names(cutv) <- colnames(data)
  #send results to outfile
  write.csv(cutv, paste("QCutoff_", name, ".csv",sep=""),row.names=F)
}

###############################################
createDAFSfile <- function(dataset, prefix) {
  ##preprocessing the data ## only need to be done once
  data <- dataset
  out <- apply(data,1,function(x) all(x==0))
  if(length(out[out=="TRUE"])>0) data <- data[-which(out=="TRUE"),]

  qcut(data, prefix) ### OBS 
  co <- read.csv("QCutoff_test.csv", sep="\t")
  co_un <- 2^co  # return back to normal numbers
  
  thr <- matrix(nrow=dim(dataset)[1], ncol=dim(dataset)[2])
  # check per gene per sample if passes
  for(i in 1:dim(co)[1]){
    thr[,i] <-dataset[,i] > co_un[i]
  }
  
  keep <- rowSums(thr)/max(rowSums(thr))*100 >= 10 ## filter out genes that show an expression of less then 10%
  dataset <- dataset[keep,]
  write.table(dataset,file=paste("Results/filtered_counts_", prefix, ".txt", sep=""),quote=F,col.names=NA,row.names=T,sep="\t")
}

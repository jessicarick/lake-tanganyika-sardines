# general packages
library(tidyverse)
library(data.table)
library(reshape2)

# graphics packages
library(ggpubr)
library(ggsci)
library(viridis)

# genetics packages
library(vcfR)
library(genetics)
library(dartR)
library(adegenet)

### ----- function to calc pca from genetic covariances --------- 
do.pca<-function(gmat){  ### inds in columns, loci in rows
  gmn<-apply(gmat,1,mean, na.rm=T)  
  gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat))
  gprime<-gmat-gmnmat ## remove mean
  gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i==j){
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
      }
      else{
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        gcovarmat[j,i]<-gcovarmat[i,j]
      }
    }
  }
  prcomp(~.,data=as.data.frame(gcovarmat),center=TRUE,scale=FALSE,na.action=na.exclude)
}



### ------------------------- 

## custom function to randomize genotypes to create null distribution
## gl.obj needs to be a genlight object
## pop is a vector of population IDs for the individuals in gl.obj
## niter is the number of randomizations to run

randomize.dapc <- function(gl.obj,pop,npca,nda=3,niter=10,verbose=TRUE,return.all=FALSE){
  gl.obj$pop <- pop
  orig.mat <- as.matrix(gl.obj)
  rand.loadings <- c()
  
  for (i in 1:niter){
    if(verbose==TRUE){
      print(paste("running iteration ",i,sep=""))
    }
    rand.mat <- apply(orig.mat, 2, function(x) sample(x, replace = T))
    rand.gl <- as.genlight(rand.mat)
    rm(rand.mat)
    indNames(rand.gl) <- indNames(gl.obj)
    ploidy(rand.gl) <- 2
    rand.gl$pop <- gl.obj$pop
    
    # remove any NA loci
    toRemove <- is.na(glMean(rand.gl, alleleAsUnit = FALSE)) # TRUE where NA
    which(toRemove) # position of entirely non-typed loci
    glNoNA <- rand.gl[, !toRemove]
    
    dapc <- dapc(glNoNA,pop=pop(glNoNA),
                 n.pca=npca, n.da=nda)
    loadings <- dapc$var.contr[,1]
    rm(rand.gl)
    rm(dapc)
    #hist(loadings)
    #print(summary(loadings))
    
    rand.loadings <- append(rand.loadings,loadings)
    
    #    if(return.all == TRUE){
    #      return(rand.loadings)
    #    }
  }
  rand.hist <- hist(rand.loadings)
  #quantile(rand.loadings,c(0.95,0.975,0.99,1),na.rm=TRUE)
  if(return.all == TRUE){
    return(rand.loadings)
    return(rand.hist)
  } else {
    return(quantile(rand.loadings,c(0.99),na.rm=TRUE))
    return(rand.hist)
  }
}


### ------------------------- 

## Reich FST estimator (Reich et al. 2009; https://doi.org/10.1038/nature08365)
## vectorized version
## input=genlight object
## FST will be calculated between pops in genlight object
## specify number of bootstraps using "bootstrap=100"
## or "bootstrap=FALSE" to not calculate bootstraps

reich.fst <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    install.packages("dplyr")
    library(dplyr, character.only=T)
  }
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))
  
  fsts <- matrix(nrow=npop,
                 ncol=npop,
                 dimnames=list(levels(gl@pop),levels(gl@pop)))
  
  if (bootstrap != FALSE){
    n.bs <- bootstrap
    bs <- data.frame(matrix(nrow=nrow(combinat::combn2(levels(gl@pop))),
                            ncol=n.bs+5))
  }
  
  k <- 0
  
  for (p1 in levels(gl@pop)){
    for (p2 in levels(gl@pop)){
      if (which(levels(gl@pop) == p1) < which(levels(gl@pop) == p2)) {
        k <- 1+k
        
        pop1 <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[p2,p1] <- F
        if (verbose == TRUE) {
          print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
        }
        
        if (bootstrap != FALSE) {
          if (verbose == TRUE) {
            print("beginning bootstrapping")
          }
          
          bs[k,1:3] <- c(p2,p1,as.numeric(F))
          
          for (i in 1:n.bs){
            loci <- sample((1:nloc), nloc, replace=TRUE)
            
            pop1.bs <- matrix(as.matrix(pop1)[,loci],
                              ncol=length(loci))
            pop2.bs <- matrix(as.matrix(pop2)[,loci],
                              ncol=length(loci))
            
            a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
            a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
            n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
            n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
            
            h1 <- (a1*(n1-a1))/(n1*(n1-1))
            h2 <- (a2*(n2-a2))/(n2*(n2-1))
            
            N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
            D <- N + h1 + h2
            
            F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
            bs[k,i+5] <- F.bs
          }
          if (verbose == TRUE){
            print(paste("bootstrapping 95% CI: ",
                        quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),"-",
                        quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),
                         quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T))
        }
        
      }
    }
  }
  
  fsts[fsts < 0] <- 0
  
  if (bootstrap != FALSE){
    colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","min_CI","max_CI")
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
    
    if (plot == TRUE){
      print("drawing plot with bootstrap CI")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- bs[,1:5]
      plot.data$fst_estimate <- as.numeric(plot.data$fst_estimate)
      plot.data$min_CI <- as.numeric(plot.data$min_CI)
      plot.data$max_CI <- as.numeric(plot.data$max_CI)
      plot.data$pop_pair <- paste(plot.data$pop1,plot.data$pop2,sep="_")
      plot.data$signif <- case_when(plot.data$min_CI > 0 ~ TRUE,
                                    TRUE ~ FALSE)
      
      
      bs.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate,col=signif)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=min_CI,ymax=max_CI),width=0.1,size=1) + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(bs.plot)
    }
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
    
    if (plot == TRUE){
      print("drawing plot without bootstrap CI")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- data.frame(combinat::combn2(row.names(fsts)),
                              fst_estimate=fsts[lower.tri(fsts)])
      plot.data$pop_pair <- paste(plot.data$X1,plot.data$X2,sep="_")
      
      fst.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(fst.plot)
    }
  }
  
  return(fst.list)
  beepr::beep()
}

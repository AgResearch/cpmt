#-------------------------------------------
#two distance functions in R
#TO-DO: Julia version  
#------------------------------------------
#a wrapper to phicoef in package GenomicRanges
phicor=function(geno, return.dist=TRUE){
  require(GenomicRanges)
  nsample=nrow(geno)
  nsnp=ncol(geno)
  phicorr <- matrix(0, nrow=nsample, ncol=nsample)
  for(i in 1:nsample){
    for(j in 1:nsample){
      if(i>j)
        phicorr[i, j] = phicoef(as.logical(geno[i,]), as.logical(geno[j,]))
    }
  }
  if(return.dist==FALSE) return(phicorr)
  as.dist(1-phicorr)
}#>>>

#Hamming distance
hammingDist=function(geno){
  nsample=nrow(geno)
  nsnp=ncol(geno)
  distHamming <- matrix(0, nrow=nsample, ncol=nsample)
  for(i in 1:nsample){
    for(j in 1:nsample){
      if(i>j)
        distHamming[i, j] = sum(geno[i,] != geno[j,])
    }
  }
  as.dist(distHamming/nsnp)
}#>>>


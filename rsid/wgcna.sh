#!/usr/bin/bash

Rscript -e '
# Weighted Correlation Network Analysis
  library(dplyr)
  INF <- Sys.getenv("INF")
  outfile <- file.path(INF,"INTERVAL","o5000-inf1-outlier_in-r2.sample")
  header <- read.table(outfile, as.is=TRUE, header=TRUE, nrows=1)
  d <- read.table(outfile,skip=2,as.is=TRUE,col.names=names(header))
  prot <- d[grepl("__",names(d))]
  prot <- filter(prot,!is.na(apply(prot,1,sum)))
  names(prot) <- unlist(lapply(strsplit(names(prot),"___"),"[",1))
  pqtl <- read.table(file.path("work","INF1.merge"),header=TRUE)
  prot <- select(prot,gsub("4E","X4E",unique(pqtl[["prot"]])))
  require(WGCNA)
  enableWGCNAThreads()
# Adjacency matrix using soft thresholding with beta=6
  ADJ <- abs(cor(prot, method="pearson"))^6
# genes < 5,000
  k <- as.vector(apply(ADJ,2,sum,na.rm=TRUE))
# genes > 5,000
# k <- softConnectivity(datE=prot,power=6)
# network analysis on 70 most connected genes
# prot[,rank(-k,ties.method="first") <= 70]
# histogram of k and a scale free topology plot
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  hist(k)
  scaleFreePlot(k, main="Check scale free topology\n")
# dissimilarity Topological Overlap Matrix
  dissADJ <- 1 - ADJ
  dissTOM <- TOMdist(ADJ)
  collectGarbage()
# partition around medoids (PAM) based on dissimilarity
  require(cluster)
  pam4 <- pam(as.dist(dissADJ),4)
  pam5 <- pam(as.dist(dissADJ),5)
  pam6 <- pam(as.dist(dissADJ),6)
  pamTOM4 <- pam(as.dist(dissTOM),4)
  pamTOM5 <- pam(as.dist(dissTOM),5)
  pamTOM6 <- pam(as.dist(dissTOM),6)
  table(pam4$clustering,pamTOM4$clustering)
  table(pam5$clustering,pamTOM5$clustering)
  table(pam5$clustering,pamTOM6$clustering)
# average linkage hierachical clusterin
# ADJ
  hierADJ <- hclust(as.dist(dissADJ),method="average")
  colorStaticADJ <- as.character(cutreeStaticColor(hierADJ,cutHeight=.99,minSize=5))
  colorDynamicADJ <- labels2colors(cutreeDynamic(hierADJ,method="tree",minClusterSize=5))
  colorDynamicHybridADJ <- labels2colors(cutreeDynamic(hierADJ,distM=dissADJ,cutHeight=0.998,
                                         deepSplit=2,pamRespectsDendro=FALSE))
  colorADJ <- data.frame(pam5$clustering,colorStaticADJ,colorDynamicADJ,colorDynamicHybridADJ)
  sizeGrWindow(10,5)
  plotDendroAndColors(dendro=hierADJ,colors=colorADJ,
                      dendroLabels=FALSE,
                      marAll=c(0.2,8,2.7,0.2),
                      main="Gene dendrogram and module colors")
# TOM
  hierTOM <- hclust(as.dist(dissTOM),method="average");
  colorStaticTOM <- as.character(cutreeStaticColor(hierTOM,cutHeight=.99,minSize=5))
  colorDynamicTOM <- labels2colors(cutreeDynamic(hierTOM,method="tree",minClusterSize=5))
  colorDynamicHybridTOM <- labels2colors(cutreeDynamic(hierTOM,distM=dissTOM,cutHeight=0.998,
                                         deepSplit=2,pamRespectsDendro=FALSE))
  colorTOM <- data.frame(pamTOM5$clustering,colorStaticTOM,colorDynamicTOM,colorDynamicHybridTOM)
  sizeGrWindow(10,5)
  plotDendroAndColors(hierTOM,colors=colorTOM,
                      dendroLabels=FALSE,
                      marAll=c(1,8,3,1),
                      main="Gene dendrogram and module colors, TOM dissimilarity")
  options(width=200)
  colorADJTOM <- cbind(colorADJ,colorTOM)
  table(colorADJTOM$colorStaticADJ)
  for(col in c("blue","brown","grey","turquoise")) print(subset(colorADJTOM,colorStaticADJ==col))
'

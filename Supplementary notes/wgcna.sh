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
# Ajacency matrix using soft thresholding with beta=6
  ADJ <- abs(cor(prot, method="pearson"))^6
# genes < 5,000
  k <- as.vector(apply(ADJ,2,sum,na.rm=TRUE))
# genes > 5,000
# k <- softConnectivity(datE=out,power=6)
# histogram of k and a scale free topology plot
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  hist(k)
  scaleFreePlot(k, main="Check scale free topology\n")
# network analysis on 70 most connected genes
  datExpr=prot[,rank(-k,ties.method="first") <= 70]
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
  hierADJ <- hclust(as.dist(dissADJ),method="average")
# dendrogram
  sizeGrWindow(10,5)
  plotDendroAndColors(hierADJ,colors=pam4$clustering,dendroLabels=FALSE,hang=0.03)
  title("Gene hierarchical clustering dendrogram")
  colorStaticADJ <- as.character(cutreeStaticColor(hierADJ,cutHeight=.99,minSize=10))
# dendrogram with module colors
  sizeGrWindow(10,5)
  plotDendroAndColors(hierADJ,colors=data.frame(pam4$clustering,colorStaticADJ),dendroLabels=FALSE,abHeight=0.99,
                      main="Gene dendrogram and module colors")
# module definition via dynamic branch cutting methods
  branch.number <- cutreeDynamic(hierADJ,method="tree")
# transforms the branch numbers into colors
  colorDynamicADJ <- labels2colors(branch.number)
  colorDynamicHybridADJ <- labels2colors(cutreeDynamic(hierADJ,distM=dissADJ,cutHeight=0.998,deepSplit=2,pamRespectsDendro=FALSE))
# all module detection methods together
  sizeGrWindow(10,5)
  plotDendroAndColors(dendro=hierADJ,colors=data.frame(pam4$clustering,colorStaticADJ,colorDynamicADJ,colorDynamicHybridADJ),
                      dendroLabels=FALSE,marAll=c(0.2,8,2.7,0.2),main="Gene dendrogram and module colors")
# dendrogram
  hierTOM <- hclust(as.dist(dissTOM),method="average");
  colorStaticTOM <- as.character(cutreeStaticColor(hierTOM,cutHeight=.99,minSize=20))
  colorDynamicTOM <- labels2colors(cutreeDynamic(hierTOM,method="tree"))
  colorDynamicHybridTOM <- labels2colors(cutreeDynamic(hierTOM,distM=dissTOM,cutHeight=0.998,
                                         deepSplit=2,pamRespectsDendro = FALSE))
  sizeGrWindow(10,5)
  plotDendroAndColors(hierTOM,colors=data.frame(hierTOM$labels,colorStaticTOM,colorDynamicTOM,colorDynamicHybridTOM), 
                      dendroLabels=FALSE, marAll=c(1,8,3,1),
                      main="Gene dendrogram and module colors, TOM dissimilarity")
'

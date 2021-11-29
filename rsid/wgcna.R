#!/usr/bin/bash

Rscript -e '
# Data handling
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
# Weighted Correlation Network Analysis
  suppressMessages(require(WGCNA))
  enableWGCNAThreads()
# Adjacency matrix using soft thresholding with beta=6
  ADJ <- abs(cor(prot, method="pearson"))^6
# histogram of k and a scale free topology plot
  k <- as.vector(apply(ADJ,2,sum,na.rm=TRUE))
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
  for(j in 4:6)
  {
    pam_name <- paste0("pam",j)
    pamTOM_name <- paste0("pamTOM",j)
    assign(pam_name, pam(as.dist(dissADJ),j))
    assign(pamTOM_name,pam(as.dist(dissTOM),j))
    tc <- table(get(pam_name)$clustering,get(pamTOM_name)$clustering)
    print(tc)
    print(diag(tc))
  }
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
# Further correlations
  corRaw <- cor(prot)
  diag(corRaw) <- 0
  distance <- as.dist(1-abs(corRaw))
  colnames(corRaw) <- rownames(corRaw) <- names(prot)
  suppressMessages(require(reshape))
  r <- melt(corRaw) %>% mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  colorADJTOM_nogrey <- subset(colorADJTOM,colorStaticTOM!="grey")
  r_nogrey <- melt(corRaw[rownames(colorADJTOM_nogrey),rownames(colorADJTOM_nogrey)]) %>% mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  library(RCy3)
  cytoscapePing()
  cytoscapeVersionInfo ()
  nodes <- data.frame(id=gsub("X4","4",rownames(colorADJTOM_nogrey)),
           group=with(colorADJTOM_nogrey,colorStaticTOM),
           stringsAsFactors=FALSE)
  edges <- data.frame(source=with(r_nogrey,gsub("X4","4",X1)),
           target=with(r_nogrey,gsub("X4","4",X2)),
           weight=with(r_nogrey,value),
           stringsAsFactors=FALSE) %>% filter(!is.na(weight))
  createNetworkFromDataFrames(nodes,edges,title="protein network", collection="DataFrame")
  nodedata <- getTableColumns("node")
  selectNodes(subset(nodedata,group=="turquoise")$name, by='id', pre=FALSE)
  createSubnetwork(subset(nodedata,group=="turquoise")$name,"name")
  exportImage("turquoise.png",type="PNG",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
  saveSession("turquoise.cys")
  suppressMessages(require(Biobase))
  suppressMessages(library(GOstats))
  gData <- new("ExpressionSet", exprs=t(prot))
  corrGraph = compCorrGraph(gData, k=6, tau=0.6)
  edgemode(corrGraph) <- "undirected"
  plot(corrGraph)
  createNetworkFromGraph(corrGraph,"myGraph")
  exportImage("corrGraph.png",type="PNG",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
  require(igraph)
  g <- graph_from_graphnel (corrGraph)
  plot(g)
  write_graph(g,"igraph.el","edgelist")
# PW-pipeline codes
  library(diagram)
  sel <- with(nodedata,name)
  colnames(corRaw) <- rownames(corRaw) <- gsub("X4E","4E",colnames(corRaw))
  plotmat(round(corRaw[sel,sel],2))
  require(network)
  n <- network(m, directed=FALSE)
  plot(n)
  require(graph)
  gmat <- new("graphAM", adjMat=m, edgemode='undirected')
  glist <- as(gmat, 'graphNEL')
  plot(glist)
  require(Rgraphviz)
# genes > 5,000
# k <- softConnectivity(prot,power=6)
# network analysis on 70 most connected genes
# prot[,rank(-k,ties.method="first") <= 70]
  library(Rtsne)
  rtsne <- Rtsne(as.matrix(prot),dims=3,perplexity=15,theta=0.25,pc=FALSE)
  plot(rtsne$Y,asp=1)
'

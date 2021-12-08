#!/rds/user/jhz22/hpc-work/bin/Rscript --vanilla

library(dplyr)
library(pQTLtools)
target <- inf1["target.short"]
rownames(target) <- inf1[["prot"]]
INF <- Sys.getenv("INF")
INF_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
             left_join(inf1,by="prot")
rsid <- unique(INF_METAL[["rsid"]])
rsid_prot <- with(INF_METAL,cbind(rsid,target.short))
cis <- subset(INF_METAL,cis.trans=="cis")
trans <- subset(INF_METAL,cis.trans=="trans")
outfile <- file.path(INF,"INTERVAL","o5000-inf1-outlier_in-r2.sample")
header <- read.table(outfile, as.is=TRUE, header=TRUE, nrows=1)
d <- read.table(outfile,skip=2,as.is=TRUE,col.names=names(header))
prot <- d[grepl("__",names(d))]
prot <- filter(prot,!is.na(apply(prot,1,sum)))
names(prot) <- unlist(lapply(strsplit(gsub("X4E","4E",names(prot)),"___"),"[",1))
prot <- select(prot,INF_METAL[["prot"]])
names(prot) <- target[names(prot),1]

library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()
deleteAllNetworks()
suppressMessages(require(Biobase))
suppressMessages(library(GOstats))
gData <- new("ExpressionSet", exprs=t(prot))
corrGraph = compCorrGraph(gData, k=6, tau=0.7)
edgemode(corrGraph) <- "undirected"
plot(corrGraph)
suid_corrGraph <- createNetworkFromGraph(corrGraph,"corrGraph")
addCyNodes(rsid)
sapply(1:nrow(rsid_prot),function(x) addCyEdges(rsid_prot[x,]))
layoutNetwork("attribute-circle")
exportImage(file.path(INF,"Cytoscape","corrGraph.pdf"),type="PDF",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
exportNetwork(file.path(INF,"Cytoscape","corrGraph.sif"))
saveSession(file.path(INF,"Cytoscape","corrGraph.cys"),overwriteFile=TRUE)
deleteNetwork(suid_corrGraph)

require(igraph)
g <- graph_from_graphnel(corrGraph) +
     vertices(unique(cis[["rsid"]]),color="red") +
     vertices(unique(trans[["rsid"]]),color="blue") + edges(as.vector(t(rsid_prot)))
plot(g)
suid_corrIGraph <- createNetworkFromIgraph(g,"corrIGraph")
layoutNetwork("attribute-circle")
exportImage(file.path(INF,"Cytoscape","corrIGraph.pdf"),type="PDF",resolution=300,height=8,width=12,units="in",overwriteFile=TRUE)
exportNetwork(file.path(INF,"Cytoscape","corrIGraph.sif"))
saveSession(file.path(INF,"Cytoscape","corrIGraph.cys"),overwriteFile=TRUE)
deleteNetwork(suid_corrIGraph)

wgcna_etc <- function()
# Weighted Correlation Network Analysis
{
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
  colorTOM <- data.frame(pamTOM5$clustering,colorStaticTOM,colorDynamicTOM)
  pdf(file.path(INF,"Cytoscape","pamTOM.pdf"))
  plotDendroAndColors(hierTOM,colors=colorTOM,
                      dendroLabels=FALSE,
                      marAll=c(1,8,3,1),
                      main="Gene dendrogram and module colors, TOM dissimilarity")
  dev.off()
  options(width=200)
  colorADJTOM <- cbind(colorADJ,colorTOM)
  table(colorADJTOM$pamTOM5.clustering)
  for(x in 1:5) print(subset(colorADJTOM,pamTOM5.clustering==x))
  table(colorADJTOM$colorDynamicTOM)
  for(col in c("blue","brown","grey","turquoise","yellow")) print(subset(colorADJTOM,colorDynamicTOM==col))
# Further correlations
  corRaw <- cor(prot)
  diag(corRaw) <- 0
  distance <- as.dist(1-abs(corRaw))
  colnames(corRaw) <- rownames(corRaw) <- names(prot)
  suppressMessages(require(reshape))
  r <- melt(corRaw) %>% mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  colorADJTOM_nogrey <- subset(colorADJTOM,colorStaticTOM!="grey")
  r_nogrey <- melt(corRaw[rownames(colorADJTOM_nogrey),rownames(colorADJTOM_nogrey)]) %>%
              mutate(value=ifelse(X1!=X2 & value>=0.7,value,NA))
  nodes <- data.frame(id=gsub("X4","4",rownames(colorADJTOM_nogrey)),
           group=with(colorADJTOM_nogrey,colorStaticTOM),
           stringsAsFactors=FALSE)
  edges <- data.frame(source=with(r_nogrey,gsub("X4","4",X1)),
           target=with(r_nogrey,gsub("X4","4",X2)),
           weight=with(r_nogrey,value),
           stringsAsFactors=FALSE) %>% filter(!is.na(weight))
  suid_wgnca <- createNetworkFromDataFrames(nodes,edges,title="turquoise", collection="DataFrame")
  layoutNetwork("attribute-circle")
  nodedata <- getTableColumns("node")
  selectNodes(subset(nodedata,group=="turquoise")$name, by='id', pre=FALSE)
  createSubnetwork(subset(nodedata,group=="turquoise")$name,"name")
  exportImage(file.path(INF,"Cytoscape","turquoise.pdf"),type="PDF",overwriteFile=TRUE)
  exportNetwork(file.path(INF,"Cytoscape","turquoise.sif"))
  saveSession(file.path(INF,"Cytoscape","turquoise.cys"))
}

wgcna_etc()

library(RColorBrewer)
string.cmd = 'string disease query disease="multiple sclerosis" cutoff=0.9 species="Homo sapiens" limit=10000'
commandsRun(string.cmd)
getTableColumnNames('node')
Nodes <- getAllNodes()
ENSP <- data.frame(ensp=gsub("9606.","",Nodes))
ENS <- read.table(file.path(INF,"work","ensGtp.txt.gz"),col.names=c("ensg","enst","ensp"),sep="\t")
ENST <- read.table(file.path(INF,"work","ensemblToGeneName.txt.gz"),col.names=c("enst","symbol"))

ms <- function(d,filename)
{
  print(dim(d))
  inf1_nodes <- with(d,paste0("9606.",ensp))
  Names <- getTableColumns('node',"name") %>% filter(name %in% inf1_nodes) %>% rownames()
  suid_INF1 <- createSubnetwork(Names,subnetwork.name=filename)
  suid_INF1connected <- createSubnetwork(edges='all',subnetwork.name=paste0(filename,"connected"))
  layoutNetwork("attribute-circle")
  exportNetwork(file.path(INF,"Cytoscape",paste0(filename,".sif")))
  saveSession(file.path(INF,"Cytoscape",paste0(filename,".cys")))
}
d <- left_join(ENSP,ENS) %>% left_join(ENST) %>% left_join(inf1, by=c('ensg'='ensembl_gene_id')) %>% filter(symbol==gene)
ms(d,"MS")
d <- filter(d, prot %in% INF_METAL$prot)
ms(d,"MS-70")
library(ndexr)
ndexcon <- ndex_connect()
networks <- ndex_find_networks(ndexcon, "Multiple Sclerosis")
print(networks[,c("name","nodeCount","edgeCount")])
networks <- ndex_find_networks(ndexcon, "Trastuzumab")
print(networks[,c("name","nodeCount","edgeCount")])
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)
suid_trastuzumab <- importNetworkFromNDEx(networkId)

# ensGtp.txt.gz
# ENSG00000215700	ENST00000400776	ENSP00000383587
# ensemblToGeneName.txt.gz
# ENST00000400776	PNRC2

misc <- function()
  library(diagram)
  sel <- with(nodedata,name)
  colnames(corRaw) <- rownames(corRaw) <- gsub("X4E","4E",colnames(corRaw))
  plotmat(round(corRaw[sel,sel],2))
  require(network)
  m <- abs(corRaw-diag(corRaw))
  n <- network(m, directed=FALSE)
  plot(n)
  require(graph)
  gmat <- new("graphAM", adjMat=m, edgemode='undirected')
  glist <- as(gmat, 'graphNEL')
  plot(glist)
# genes > 5,000
# k <- softConnectivity(prot,power=6)
# network analysis on 70 most connected genes
# prot[,rank(-k,ties.method="first") <= 70]
  library(Rtsne)
  rtsne <- Rtsne(as.matrix(prot),dims=3,perplexity=15,theta=0.25,pca=FALSE)
  plot(rtsne$Y,asp=1)
}

obsolete <- function()
{
  pheatmap(rxc,
           color = col,
           legend = TRUE,
           main = "",
           angle_col = "315",
           filename = "INF1_pQTL_immune_qtl.png",
           width = 17,
           height = 11,
           treeheight_row = 100,
           treeheigh_col = 100,
           cellheight = 20,
           cellwidth = 20,
           fontsize_row = 14,
           fontsize = 13)

  library(highcharter)
  fntltp <- JS("function(){
    return this.series.xAxis.categories[this.point.x] + ' ' +
    this.series.yAxis.categories[this.point.y] + ':<br>' +
    Highcharts.numberFormat(this.point.value, 2);
  }")
  hc <- data.frame()
  i <- 1
  for(cn in colnames(tbl)) for(rn in rownames(tbl)) {
     hc[i,c("f1","f2","v")] <- c(cn,rn,tbl[rn,cn])
     i <- i + 1
  }
  n <- 4
  stops <- data.frame(
    q = 0:n/n,
    c = c("#4287f5","grey","#ffffff","grey","#e32222"),
    stringsAsFactors = FALSE
  )

  hc$f1 <- as.factor(hc$f1)
  hc$f2 <- as.factor(hc$f2)
  f1 <- levels(hc$f1)
  highchart() %>%
    hc_title(text = "pQTLs and Immune-related Diseases",align="center")%>%
    hc_xAxis(categories = f1) %>%
    hc_yAxis(categories = hc$f2, reversed = TRUE)%>%
    hc_colorAxis(min = -1, max=1, stops=list_parse2(stops)) %>%
    hc_legend(align = "right",layout = "vertical",
              margin = 0,verticalAlign = "top",
              y = 30,symbolHeight = 200) %>%
    hc_tooltip(formatter = fntltp) %>%
    hc_add_series(data = hc, type = "heatmap",
                  hcaes(x = f1,y = f2,value = v),
                  dataLabels = list(enabled = FALSE))

  library(gap)
  aux <- with(with(mat, cbind(inv_chr_pos_a1_a2(MarkerName)[c("chr","pos")],rsid,Allele1,Allele2,prots,HLA,cistrans,efoTraits,qtl_direction)), {
              flag <- (HLA==1)
  # a bit too specific here as they involve too many proteins nevertheless each combination with the same effect allele
              Allele1[8:13] <- "T"
              Allele2[8:13] <- "C"
              Allele1[25] <- "C"
              Allele2[25] <- "G"
              colId <- paste0(substr(chr,4,5),":",pos,"(",Allele1,"/",Allele2,")")
              colId[flag] <- paste0(colId[flag],"*")
              colLabel <- paste0(colId," (",prots,")")
              col <- rep("blue",nrow(mat))
              col[cistrans=="cis"] <- "red"
              data.frame(colLabel,col,efoTraits,qtl_direction)
            })
  Col <- unique(aux[c("colLabel","col")])
  rownames(Col) <- with(Col,colLabel)

  RXC <- with(aux,table(efoTraits,colLabel))
  indices <- aux[c("efoTraits","colLabel","qtl_direction")]
  for(cn in colnames(RXC)) for(rn in rownames(RXC)) {
     s <- subset(indices,efoTraits==rn & colLabel==cn)
     qd <- s[["qtl_direction"]]
     if(length(qd)>1) stop("duplicates")
     class(qd) <- "numeric"
     if(nrow(s)>0 & !is.na(qd[1])) RXC[rn,cn] <- qd[1]
  }

  library(gplots)
  png("INF1_pQTL_immune_gplots.png",height=35,width=40,units="cm",res=300)
  heatmap.2(RXC, scale = "none", keysize=0.8, col = colorpanel(5, "blue", "white", "red"), margin=c(20,20), trace = "none",
            colCol=Col[colnames(RXC),"col"], dendrogram="none", density.info = "none", srtCol=45)
  dev.off()

  efo_list_immune <- subset(read.csv("work/efo_list_annotated.csv",as.is=TRUE),immune_mediated==1)
  isd1 <- merge(aggr,subset(ps,efo%in%with(efo_list_immune,EFO)),by="hg19_coordinates")
  write.table(isd1,file="isd1.tsv",row.names=FALSE,quote=FALSE,sep="\t")

  load(file.path(INF,"files","efo.rda"))
  efo_0000540 <- gsub(":","_",as.data.frame(isd)[["efo_0000540"]])
  isd2 <- merge(aggr,subset(ps,efo%in%efo_0000540),by="hg19_coordinates")
  write.table(isd2,file="isd2.tsv",row.names=FALSE,quote=FALSE,sep="\t")

  fang_efo <- gsub(":","_",with(read.delim("doc/fang.efos.txt",as.is=TRUE),id))
  isd3 <- merge(aggr,subset(ps,efo%in%fang_efo),by="hg19_coordinates")
  write.table(isd3,file="isd3.tsv",row.names=FALSE,quote=FALSE,sep="\t")

# A test of colorRampPalette
  YlOrBr <- c("#4287f5","grey","#ffffff","grey","#e32222")
  filled.contour(volcano,color.palette = colorRampPalette(YlOrBr, space = "Lab"), asp = 1)
# Colouring for the dendrogram
  library(dendextend)
  Rowv <- rxc %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
     set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
     ladderize
  Colv <- rxc %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
     set("branches_k_color", k = 2, value = c("orange", "blue")) %>%
     set("branches_lwd", 1.2) %>%
     ladderize
# stats
  heatmap(scale(rxc), scale = "none")
  heatmap(scale(rxc), Rowv = Rowv, Colv = Colv, scale = "none")
# gplots
  heatmap.2(scale(rxc), scale = "none", col = bluered(100), Rowv = Rowv, Colv = Colv, trace = "none", density.info = "none")
  options(width=200)
  IL_12Bcis <- ieugwasr::phewas("rs10076557")
  data.frame(IL_12Bcis)
}

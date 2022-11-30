  export rt=${INF}/mr/gsmr
  (
    cat ${rt}/mr/*mr | head -1
    ls ${rt}/mr/*mr | grep -v TNFB | xargs -l -I {} sed '1d' {}
  ) > $rt/mr-efo.mr
  (
    cat ${rt}/mr/*het | head -1
    ls ${rt}/mr/*het | grep -v TNFB | xargs -l -I {} sed '1d' {}
  ) > $rt/mr-efo.het
  Rscript -e '
    library(dplyr)
    library(openxlsx)
    rt <- Sys.getenv("rt")
    IVW <- read.delim(file.path(rt,"mr-efo.mr")) %>%
           mutate(fdr=p.adjust(pval,method="fdr")) %>%
           left_join(select(pQTLtools::inf1,prot,gene)) %>%
           select(gene,id.outcome,b,se,pval,fdr,nsnp,disease)
    cat(nrow(filter(IVW,fdr<=0.05)),"\n")
    write.table(IVW,file.path(rt,"mr-efo-mr.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
    Heterogeneity <- read.delim(file.path(rt,"mr-efo.het")) %>%
                     mutate(fdr=p.adjust(Q_pval,method="fdr")) %>%
                     left_join(select(pQTLtools::inf1,prot,gene)) %>%
                     select(gene,id.outcome,method,Q,Q_df,Q_pval,fdr,disease)
    cat(nrow(filter(Heterogeneity,fdr<=0.05)),"\n")
    write.table(Heterogeneity,file.path(rt,"mr-efo-het.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
#   GSMR <- read.delim(file.path(rt,"out","gsmr-efo.txt")) %>%
#           rename(target.short=protein) %>%
#           left_join(select(pQTLtools::inf1,target.short,gene)) %>%
#           select(gene,MRBASEID,trait,bxy,se,p,nsnp,fdr,Ncases,Ncontrols,id,uri,Zhengetal)
    GSMR <- read.delim(file.path(rt,"gsmr-efo-reduce.txt")) %>%
            rename(target.short=protein) %>%
            left_join(select(pQTLtools::inf1,target.short,gene))
    xlsx <- file.path(rt,"gsmr-mr.xlsx")
    wb <- createWorkbook(xlsx)
    hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
    for (sheet in c("GSMR","IVW","Heterogeneity"))
    {
      addWorksheet(wb,sheet,zoom=150)
      writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
                fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
      body <- get(sheet)
      writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
      freezePane(wb, sheet, firstCol=TRUE, firstActiveRow=3)
      width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
    # width_vec_header <- nchar(colnames(body))+2
      setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
      writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
    }
    saveWorkbook(wb, file=xlsx, overwrite=TRUE)
  '

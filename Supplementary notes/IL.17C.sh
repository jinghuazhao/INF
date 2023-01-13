#!/usr/bin/bash

# gunzip -c ${INF}/METAL/IL.17C-1.tbl.gz | cut -f1-3,10-11 | gzip -f > ${INF}/work/INF-IL.17C.gz

Rscript -e '
  suppressMessages(library(dplyr))
  suppressMessages(library(gap))
  INF <- Sys.getenv("INF")
  genes <- data.frame(MarkerName=c("chr16:88684495_G_T"),gene=c("IL17C"),color=c("red"))
  IL.17C <- read.delim(file.path(INF,"work","INF-IL.17C.gz"),as.is=TRUE) %>%
            mutate(Z=Effect/StdErr,P=as.numeric(pvalue(Z)),color=NA) %>%
            filter(!is.na(Z)) %>%
            select(Chromosome,Position,MarkerName,Z,P,color) %>%
            left_join(genes) %>%
            arrange(Chromosome,Position)
  mhtdata <- select(IL.17C,Chromosome,Position,P,gene,color)
  cis <- with(mhtdata,Chromosome==16 & Position>=88704999-1e6 & Position<=88706881+1e6)
  mhtdata[cis,"color"] <- "red"
  subset(mhtdata,!is.na(gene))
  png(file.path(INF,"work","INF-IL.17C-mhtplot.png"), res=300, units="in", width=12, height=8)
    par(cex=0.8, mar=c(6,6,3,1),xpd=TRUE)
    ops <- mht.control(colors=rep(c("blue4","skyblue"),11),srt=0,yline=2.5,xline=2)
    hops <- hmht.control(data=filter(mhtdata,!is.na(gene)))
    mhtplot2(mhtdata,ops,hops,xlab="",ylab="",srt=0, cex.axis=2)
    axis(2,at=0:15)
    abline(h=-log10(5e-10),col="red")
  dev.off()
  IL.17C <- mutate(IL.17C,MarkerName=ifelse(!is.na(gene),gene,MarkerName))
  png(file.path(INF,"work","INF-IL.17C-mhtplot.trunc.png"), res=300, units="in", width=12, height=8)
    par(oma=c(0,0,0,0), mar=c(5,6.5,1,1),xpd=TRUE)
    mhtplot.trunc(subset(IL.17C,select=-color), chr="Chromosome", bp="Position", z="Z", snp="MarkerName",
                  suggestiveline=NULL, genomewideline=-log10(5e-10),
                  cex.mtext=1.2, cex.text=0.7,
                  annotatelog10P=-log10(1.16e-14), annotateTop = FALSE, highlight=with(genes,gene),
                  mtext.line=3, y.brk1=0.0001, y.brk2=0.0002, trunc.yaxis=FALSE, delta=0.1, cex.axis=1.2,
                  cex=0.5, font=2, font.axis=1, y.ax.space=1,
                  col = c("blue4", "skyblue"))
  dev.off()
'

suppressMessages(library(dplyr))
suppressMessages(library(gap))
INF <- Sys.getenv("INF")
genes <- data.frame(MarkerName=c("chr16:88684495_G_T"),gene=c("IL17C"),color=c("red"))
IL.17C <- read.delim(file.path(INF,"work","INF-IL.17C.gz"),as.is=TRUE) %>%
          mutate(Z=Effect/StdErr,P=as.numeric(pvalue(Z)),gene=NA,color=NA) %>%
          filter(!is.na(Z)) %>%
          select(Chromosome,Position,MarkerName,Z,P,color) %>%
          arrange(Chromosome,Position)
loc <- with(IL.17C,MarkerName=="chr16:88684495_G_T")
IL.17C[loc,"gene"] <- "IL17C"
IL.17C <- mutate(IL.17C,MarkerName=ifelse(!is.na(gene),gene,MarkerName))
png(file.path(INF,"work","SF-IL-17C.png"), res=300, units="in", width=12, height=8)
    par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
    mhtplot.trunc(subset(IL.17C,select=-color), chr="Chromosome", bp="Position", z="Z", snp="MarkerName",
                  suggestiveline=NULL, genomewideline=-log10(5e-10),
                  cex.mtext=1.2, cex.text=0.7,
                  annotatelog10P=-log10(1.16e-14), annotateTop = FALSE, highlight="IL17C",
                  mtext.line=3, trunc.yaxis=FALSE, delta=0.1, cex.axis=1.2,
                  cex=0.5, font=2, font.axis=1, y.ax.space=1,
                  col = c("blue4", "skyblue"))
dev.off()

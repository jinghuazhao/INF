#!/usr/bin/bash

Rscript -e '
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
  gz <- gzfile(file.path(INF,"METAL","IL.12B-1.tbl.gz"))
  IL.12B <- within(read.delim(gz,as.is=TRUE), {Z <- Effect/StdErr; P <- pvalue(Z); log10P <- -log10p(Z)}) %>%
            select(Chromosome,Position,MarkerName,Z,P,log10P)
  genes <- data.frame(chr=c("chr3","chr3","chr5","chr6","chr12","chr13","chr14","chr14"),
                      snpid=c("chr3:5026008_A_G",
                              "chr3:188115682_A_C",
                              "chr5:158792819_C_G",
                              "chr6:31154493_A_G",
                              "chr12:111884608_C_T",
                              "chr13:28604007_C_T",
                              "chr14:68760141_C_T",
                              "chr14:103230758_C_G"),
                      snp=c("rs11130215","rs9815073","rs10076557","rs3130510","rs3184504","rs76428106","rs12588969","rs1950897"),
INF
                      gene=c("BHLHE40","LPP","IL12B","MHC","SH2B3;TRAFD1","FLT3","RAD51B","TARF3")
           )
  IL.12B <- left_join(IL.12B,genes,by=c("MarkerName"="snpid"),keep=TRUE) %>%
             mutate(MarkerName=ifelse(!is.na(snpid),gene,MarkerName)) %>%
             select(-c(chr,snp,snpid))
  save(IL.12B, genes, file=file.path(INF,"work","IL.12B.rda"))
  log10p <- gap::log10p
  load(file.path(INF,"work","IL.12B.rda"))
  subset(IL.12B,!is.na(gene))
  png("IL.12B-mhtplot.trunc.png", res=300, units="in", width=9, height=6)
  par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
  source(file.path(INF,"csd3","IL.12B-mhtplot.trunc.R"))
  mhtplot.trunc(IL.12B, chr="Chromosome", bp="Position", z="Z", snp="MarkerName",
                suggestiveline=FALSE, genomewideline=-log10(5e-10),
                cex.mtext=1.2, cex.text=1.2,
                annotatelog10P=-log10(5e-10), annotateTop = FALSE, highlight=with(genes,gene),
                mtext.line=3, y.brk1=115, y.brk2=300, delta=0.01, cex.axis=1.5, cex=1.48, font=3, font.axis=1.5,
                y.ax.space=20,
                col = c("blue4", "skyblue")
  )
  dev.off()
'

# IL.12B[which(IL.12B$MarkerName%in%genes$snpid),"MarkerName"] <- genes[match(interaction(IL.12B[which(IL.12B$MarkerName%in%genes$snpid),c("MarkerName")]), 
#                                              interaction(genes[,c("snpid")])), "gene"]
# left_join(IL.12B,genes[,c(1:3)],by=c("MarkerName"="snpid"))%>%
# mutate(MarkerName=ifelse(snpid!="",gene,MarkerName))%>%
# select(-snp)
# chr3  chr3:5026008_A_G    rs11130215 BHLHE40
# chr3  chr3:188115682_A_C  rs9815073  LPP
# chr5  chr5:158792819_C_G  rs10076557 IL12B
# chr6  chr6:31154493_A_G   rs3130510  MHC
# chr12 chr12:111884608_C_T rs3184504  SH2B3;TRAFD1
# chr13 chr13:28604007_C_T  rs76428106 FLT3
# chr14 chr14:68760141_C_T  rs12588969 RAD51B
# chr14 chr14:103230758_C_G rs1950897  TARF3

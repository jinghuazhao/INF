#!/usr/bin/bash

module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL

(
cat STABILITY.list | \
parallel -j8 -C' ' '
  export protein={};
  echo {}
  R --no-save -q <<\ \ END
    protein <- Sys.getenv("protein");
    print(protein);
    gz <- gzfile(paste0("sumstats/STABILITY/STABILITY.",protein,".gz"));
    .libPaths("/services/tools/R/3.5.0/lib64/R/library:/data/jinhua/R")
    require(qqman);
    tbl <- read.delim(gz,as.is=TRUE);
    tbl <- within(tbl,{
       SNP <- SNPID
       CHR <- as.numeric(CHR)
       BP <- as.numeric(POS)
       P <- as.numeric(PVAL)
     })
     tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
   # manhattan <- paste0("STABILITY.",protein,".manhattan.png");
   # png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
   # manhattan(tbl,main=protein,genomewideline=-log10(5e-10),cex=0.8, col=c("blue","orange"),suggestiveline=FALSE,ylim=c(0,25));
   # dev.off();
     library(gap, lib.loc="/data/jinhua/R")
     cat(protein,"GC.lambda=",gc.lambda(with(tbl,P)),"\n")
  END
'
) > STABILITY.lambda

# ls sumstats/STABILITY/*gz | \
# sed 's|sumstats/STABILITY/STABILITY.||g;s/.gz//g' > STABILITY.list
# ::: IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN IL.10RA IL.5 TNF

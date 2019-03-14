# 14-3-2019 JHZ

protein <- Sys.getenv("protein");
print(protein);
gz <- gzfile(paste0("sumstats/STABILITY/STABILITY.",protein,".gz"));
.libPaths("/services/tools/R/3.5.0/lib64/R/library")
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
library(gap)
cat(protein,"GC.lambda=",gc.lambda(with(tbl,PVAL)),"\n")

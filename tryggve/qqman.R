# 10-7-2019 JHZ

protein <- Sys.getenv("protein");
print(protein);
gz <- gzfile(paste0("METAL/",protein,"-1.tbl.gz"));
R_LIBS <- Sys.getenv("R_LIBS")
.libPaths(R_LIBS)
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
tbl <- within(tbl,{
   SNP <- MarkerName
   CHR <- as.numeric(Chromosome)
   BP <- Position
   P <- 10^log.P.
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
qq <- paste0("METAL/",protein,".qq.png");
png(qq,width=12,height=10,units="in",pointsize=4,res=300)
qq(with(tbl,P))
dev.off()
manhattan <- paste0("METAL/",protein,".manhattan.png");
png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
manhattan(tbl,main=protein,genomewideline=-log10(5e-10),suggestiveline=FALSE,ylim=c(0,25));
dev.off();

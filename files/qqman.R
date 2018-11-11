# 11-11-2018 JHZ

protein <- Sys.getenv("protein");
print(protein);
gz <- gzfile(paste0("METAL/",protein,"-1.tbl.gz"));
qqman <- paste0("METAL/",protein,"-qqman.pdf");
.libPaths("/services/tools/R/3.5.0/lib64/R/library")
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
tbl <- within(tbl,{
   SNP <- MarkerName
   CHR <- as.numeric(Chromosome)
   BP <- Position
   P <- P.value
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
pdf(qqman,width=12,height=10)
qq(with(tbl,P))
manhattan(tbl,main=protein,genomewideline=-log10(5e-10),suggestiveline=FALSE,ylim=c(0,25));
dev.off();

# 10-11-2018 JHZ

protein <- Sys.getenv("protein");
print(protein);
gz <- gzfile(paste0("METAL/",protein,"-1.tbl.gz"));
qqman <- paste0(protein,"-qqman.png");
MarkerName <- "MarkerName";
PVAL <- "P.Value";
.libPaths("/services/tools/R/3.5.0/lib64/R/library")
require(qqman);
tbl <- read.delim(gz,as.is=TRUE);
chrpos_a1_a2 <- strsplit(gsub("chr","",tbl[MarkerName]),":")
tbl <- within(tbl,{
   SNP <- tbl[MarkerName]
   CHR <- as.numeric(unlist(lapply(chrpos_a1_a2,"[",1)))
   Pos_a1_a2 <- unlist(lapply(chrpos_a1_a2,"[",2))
   Pos_a1_a2 <- strsplit(Pos_a1_a2,"_")
   BP <- as.numeric(unlist(lapply(Pos_a1_a2,"[",1)))
   P <- tbl[PVAL]
})
tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
png(qqman,res=300,width=12,height=10,units="in")
par(mfrow=c(2,1))
qq(with(tbl,P))
manhattan(tbl,main=protein,genomewideline=-log10(5e-10),suggestiveline=FALSE,ylim=c(0,10));
dev.off();

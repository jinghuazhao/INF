# 29-8-2021 JHZ

require(dplyr)
INF <- Sys.getenv("INF")
protein <- Sys.getenv("protein")
file <- Sys.getenv("file")
prot <- filter(gap.datasets::inf1,prot==protein)[["target.short"]]
print(protein)
gz <- gzfile(file)
require(qqman)
tbl <- read.delim(gz,as.is=TRUE) %>%
       mutate(SNP=MarkerName, CHR=as.numeric(Chromosome), BP=Position, P=10^log.P.) %>%
       filter(!is.na(CHR)&!is.na(BP)&!is.na(P))
qqman <- file.path(INF,"plots",paste0(prot,"-qqman.png"))
png(qqman,width=12,height=10,units="in",pointsize=4,res=300)
par(mfrow=c(1,2))
qq(with(tbl,P),main=prot,cex.axis=2.5,cex.lab=2.5)
manhattan(tbl,genomewideline=-log10(5e-10),suggestiveline=FALSE,ylim=c(0,25),cex.axis=2.5,cex.lab=2.5)
dev.off();

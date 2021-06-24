options(echo=FALSE)

library(gap)
prot <- Sys.getenv("prot")
tag <- Sys.getenv("tag")
p <- read.delim(paste0("work/",prot,tag,".p"),as.is=TRUE)
chrs <- with(p,unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p,Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
 # lines 42 and 45 has r2 in gap_1.2.1 rather than as.numeric(r2)+1
  sentinels(ps, prot, 1)
}

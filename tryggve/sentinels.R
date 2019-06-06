options(echo=FALSE)
require(gap)
prot <- Sys.getenv("prot")
tag <- Sys.getenv("tag")
p <- read.delim(paste0("work/",prot,tag,".p"),as.is=TRUE)
chrs <- with(p,unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p,Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  sentinels(ps, 1)
}

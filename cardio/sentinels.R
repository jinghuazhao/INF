# 1-6-2019 JHZ

pp <- function(p,l)
{
  nr <- nrow(p)
  posmin <- pos <- p[l,"End"]
  pmin <- pvalue <- p[l,"P.value"]
  for (r in l:nr)
  {
    posmax <- pos <- p[r,"End"]
    pvalue <- p[r,"P.value"]
    if (pos <= 1e6 + posmin)
    {
      if (pvalue < pmin)  {posmax <- pos; pmin <- pvalue}
      if(r==nr) cat(posmin,"-",pos,"dist =",pos-posmin,"pmin=",pmin,"\n")
    }
    else pp(p,r)
  }
}
prot <- Sys.getenv("prot")
print(prot)
p <- read.delim(paste0(prot,".p"),as.is=TRUE)
chr <- with(p,unique(Chrom))
for(s in chr)
{
  cat("prot =",prot,"Chromosome =",s,"\n")
  ps <- subset(p,Chrom==s)
  pp(ps,1)
# print(ps[c("Chrom","End","MarkerName","P.value")])
}


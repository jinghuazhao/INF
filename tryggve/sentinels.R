# 2-6-2019 JHZ

pp <- function(p,st,debug=FALSE,flanking=1e6)
{
  nr <- nrow(p)
  p <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
  })
  if (debug) print(p[c("Chrom","End","d","s","MarkerName","P.value")])
  l <- p[1,"End"]
  u <- p[nrow(p),"End"]
  len <- with(p,max(s))
  if (len < flanking) {
    pmin <- with(p,min(P.value))
    x <- subset(p, P.value==pmin)
    m <- x[1,"End"]
    cat(prot, l, "-", u, "d =", u-l, "m =", m, "p =", pmin, "row =", row.names(x), "(case 1)\n")
  }
  else
  {
    s <- subset(p,s < flanking)
    pmin <- with(s,min(P.value))
    x <- subset(s, P.value==pmin)
    m <- x[1,"End"]
    t <- subset(p, m+flanking < End)
    qmin <- with(t,min(P.value))
    y <- subset(t, P.value==pmin)
    u <- p[nrow(t), "End"]
    if (qmin>pmin) {
      cat(prot, l, "-", u, "d =", u-l, "m =", m, "p =", pmin, "row =", row.names(x), "(case 2)\n")
      pp(p,row.names(t[nrow(t),]))
    }
    else pp(p,row.names(y[nrow(y),]))
  }
}
prot <- Sys.getenv("prot")
p <- read.delim(paste0(prot,".p"),as.is=TRUE)
chr <- with(p,unique(Chrom))
for(s in chr)
{
  cat("prot =",prot,"Chromosome =",s,"\n")
  ps <- subset(p,Chrom==s)
  pp(ps,1)
}


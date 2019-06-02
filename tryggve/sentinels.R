# 2-6-2019 JHZ

pp <- function(p,st,debug=FALSE,flanking=1e6)
{
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
  })
  if (debug) print(z[c("Chrom","End","d","s","MarkerName","P.value")])
  l <- z[1, "End"]
  u <- z[nrow(p), "End"]
  len <- with(z, max(s))
  if (len < flanking) {
    p1 <- with(z, min(P.value))
    x <- subset(z, P.value==p1)
    r1 <- row.names(x)[1]
    m <- x[1,"End"]
    cat(prot, l, "-", u, "d =", u-l, "m =", m, "p =", p1, "row =", r1, "(case 1)\n")
  }
  else
  {
    s <- subset(z, s < flanking)
    p1 <- with(s, min(P.value))
    x <- subset(s, P.value==p1)
    r1 <- row.names(x)[1]
    m <- x[1, "End"]
    t <- subset(z, End > m & End < m + flanking)
    if (nrow(t)==0) {
print(r1)
      r2 <- as.numeric(r1) + 1
      pp(p, r2)
    } else {
      p2 <- with(t, min(P.value))
      y <- subset(t, P.value==p2)
      u <- p[nrow(t), "End"]
      r2 <- row.names(t)[nrow(t)]
      if (p2 > p1) {
        cat(prot, l, "-", u, "d =", u-l, "m =", m, "p =", p1, "row =", r1, "(case 2)\n")
        if (r2 < nr) pp(p, r2)
      } else {
        r2 <- row.names(y)[nrow(y)]
        if (r2 < nr) pp(p, r2)
      }
    }
  }
}
prot <- Sys.getenv("prot")
p <- read.delim(paste0(prot,".p"),as.is=TRUE)
chr <- with(p,unique(Chrom))
for(s in chr)
{
  cat("prot =",prot,"Chromosome =",s,"\n")
  ps <- subset(p,Chrom==s)
  pp(ps, 1)
}

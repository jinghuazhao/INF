options(echo=FALSE)
pp <- function(p,st,debug=FALSE,flanking=1e+6)
{
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
  })
  if (debug) print(z[c("Chrom","End","d","s","MarkerName","P.value")])
  cat("st=",st, "nr=",nr)
  if (z[nrow(z),"s"] <= flanking) {
    l <- z[1, "End"]
    u <- z[nrow(z), "End"]
    p1 <- with(z, min(P.value))
    x <- subset(z, P.value==p1)
    r1 <- row.names(x)[1]
    m <- x[1,"End"]
    n <- x[1, "MarkerName"]
    cat(prot, chr, n, l, u, u-l, p1, r1, "I\n", sep=",")
  } else {
    s <- subset(z, s <= flanking)
    l <- s[1, "End"]
    u <- s[nrow(s), "End"]
    p1 <- with(s, min(P.value))
    x <- subset(s, P.value==p1)
    r1 <- row.names(x)[1]
    m <- x[1, "End"]
    n <- x[1, "MarkerName"]
    t <- subset(z, End > m & End < m + flanking)
    if (nrow(t)==0) {
      if (u-l > 0) cat(prot, chr, n, l, u, u-l, p1, r1, "II\n", sep=",")
      message(paste0("No variants +1 MB downstream so move to next block (",prot,")"))
      r2 <- as.numeric(r1) + 1
      pp(p, r2)
    } else {
      p2 <- with(t, min(P.value))
      y <- subset(t, P.value==p2)
      u <- p[nrow(t), "End"]
      r2 <- row.names(t)[nrow(t)]
      if (p2 > p1) {
        cat(prot, chr, n, l, u, u-l, p1, r1, "III\n", sep=",")
        if (r2 < nr) pp(p, r2)
      } else {
        r2 <- as.numeric(row.names(y)[nrow(y)])
        if(r2 < nr) pp(p, r2)
      }
    }
  }
}
prot <- Sys.getenv("prot")
p <- read.delim(paste0("work/",prot,".p"),as.is=TRUE)
chrs <- with(p,unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p,Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  pp(ps, 1)
}

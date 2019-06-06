options(echo=FALSE)

log10p <- function(z)
  log(2, base=10)+pnorm(-abs(z), lower.tail=TRUE, log.p=TRUE)/log(10)

pp <- function(p,st,debug=FALSE,flanking=1e+6)
{
  first <- TRUE
  nr <- nrow(p)
  z <- within(p[st:nr,],{
    d <- c(0,diff(End))
    s <- cumsum(d)
    log10p <- -log10p(Effect/StdErr)
  })
  if (debug) print(z[c("Chrom","End","d","s","MarkerName","P.value")])
  if (z[nrow(z),"s"] <= flanking & first) {
    l <- z[1, "End"]
    u <- z[nrow(z), "End"]
    log10p1 <- with(z, max(log10p))
    x <- subset(z, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- x[1,"End"]
    n <- x[1, "MarkerName"]
    cat(prot, n, l, u, u-l, log10p1, r1, "I\n", sep=",")
  } else {
    first <- FALSE
    s <- subset(z, s <= flanking)
    l <- s[1, "End"]
    u <- s[nrow(s), "End"]
    log10p1 <- with(s, max(log10p))
    x <- subset(s, log10p==log10p1)
    r1 <- row.names(x)[1]
    m <- x[1, "End"]
    n <- x[1, "MarkerName"]
    t <- subset(z, End > m & End < m + flanking)
    if (nrow(t)==0) {
      # cat(prot, n, l, u, u-l, log10p1, r1, "II\n", sep=",")
      message(paste0("No variants +1 MB downstream so move to next block (",prot,")"))
      r2 <- as.numeric(r1) + 1
      pp(p, r2)
    } else {
      log10p2 <- with(t, max(log10p))
      y <- subset(t, log10p==log10p2)
      u <- p[nrow(t), "End"]
      r2 <- row.names(t)[nrow(t)]
      if (log10p2 < log10p1) {
        cat(prot, n, l, u, u-l, log10p1, r1, "III\n", sep=",")
        if (r2 < nr) pp(p, r2)
      } else {
        r2 <- as.numeric(row.names(y)[nrow(y)])
        if(r2 < nr) pp(p, r2)
      }
    }
  }
}
prot <- Sys.getenv("prot")
tag <- Sys.getenv("tag")
p <- read.delim(paste0("work/",prot,tag,".p"),as.is=TRUE)
chrs <- with(p,unique(Chrom))
for(chr in chrs)
{
  ps <- subset(p,Chrom==chr)
  row.names(ps) <- 1:nrow(ps)
  pp(ps, 1)
}

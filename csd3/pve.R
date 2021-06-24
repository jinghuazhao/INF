# 14-11-2020 JHZ

require(gap)
t <- read.delim("INF1.tbl",as.is=TRUE)
tbl <- within(t, {
    prot <- sapply(strsplit(Chromosome,":"),"[",1)
    Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
})
## to obtain variance explained
tbl <- within(tbl,
{
  x2 <- (Effect/StdErr)^2
  r2 <- x2 / (N - 2 + x2)
  v <- 1 / (N - 1)
# r
# r <- sqrt(r2)
# vr <- (1 - r2)^2/ N
# Taylor expansion
# v2 <- 2 * r2^2 * (1 + 1/ (N + 1)^2)
})
s <- with(tbl, aggregate(r2,list(prot),sum))
names(s) <- c("prot", "pve")
se2 <- with(tbl, aggregate(v,list(prot),sum))
names(se2) <- c("p1","v")
m <- with(tbl, aggregate(r2,list(prot),length))
names(m) <- c("p2","m")
pve <- cbind(s,se2,m)
ord <- with(pve, order(pve))
sink("pve.dat")
print(pve[ord, c("prot","pve","v","m")], row.names=FALSE)
sink()
write.csv(tbl,file="INF1.csv",quote=FALSE,row.names=FALSE)

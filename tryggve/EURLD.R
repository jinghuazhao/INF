# 21-4-2019 JHZ

EURLD <- read.delim("tryggve/EURLD.bed")
EURLD <- within(EURLD, {
           dist <- round((End-Start)/2/1000)
           in250k <- ifelse(dist<250,0,1)
           in500k <- ifelse(dist<500,0,1)
           in10m <- ifelse(dist<10000,0,1)
         })
ord <- with(EURLD, order(dist))
with(EURLD,summary(dist))
with(EURLD,table(in250k))
with(EURLD,table(in500k))
with(EURLD,table(in10m))
EURLD[ord,]

# 11-5-2019 JHZ

EURLD <- read.delim("tryggve/EURLD.bed")
EURLD <- within(EURLD, {
           dist <- round((End-Start)/2/1000)
           within250k <- ifelse(dist<250,1,0)
           within500k <- ifelse(dist<500,1,0)
           within10m <- ifelse(dist<10000,1,0)
         })
with(EURLD,summary(dist))
with(EURLD,table(within250k))
with(EURLD,table(within500k))
with(EURLD,table(within10m))
ord <- with(EURLD, order(dist))
EURLD <- EURLD[ord,]
write.table(EURLD,file="tryggve/EURLD.tsv",quote=FALSE,row.names=FALSE,sep="\t")
EURLD

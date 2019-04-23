# 23-4-2019 JHZ

EURLD <- read.delim("tryggve/EURLD.bed")
EURLD <- within(EURLD, {
           dist <- round((End-Start)/2/1000)
           in250k <- ifelse(dist<250,1,0)
           in500k <- ifelse(dist<500,1,0)
           in10m <- ifelse(dist<10000,1,0)
         })
ord <- with(EURLD, order(dist))
with(EURLD,summary(dist))
with(EURLD,table(in250k))
with(EURLD,table(in500k))
with(EURLD,table(in10m))
write.csv(cbind(EURLD[ord,],line=1:nrow(EURLD)),file="tryggve/EURLD.csv",quote=FALSE,row.names=FALSE)

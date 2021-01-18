# 18-1-2021 JHZ

p <- Sys.getenv("p")
r <- Sys.getenv("r")
pr <- paste0(p,"-",r)

tbl <- read.delim(paste0(pr,".z"),sep=" ")
z <- gap::cs(tbl)
write.table(z[["MarkerName"]],file=paste0(pr,".cs"),col.names=FALSE,row.names=FALSE,quote=FALSE)

# tbl <- within(tbl,{logp <- logp(Effect/StdErr)})
# l <- cs(tbl,log_p="logp")

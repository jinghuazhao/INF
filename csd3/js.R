# 10-2-2020 JHZ

library(jsonlite)
SomaLogic <- read_json("js.SomaLogic")
xdata <- SomaLogic$x$data
library(gap)
d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
r <-  mhtplot2d(d)
r <- within(r,{x=x/1.3e8; y=y/1.3e8; log10p <- abs(log10p)})
cols <- c("id","chr1","x","gene","target","log10p","chr2","y","col")
d3c <- subset(r[cols],col=="blue")
d3t <- subset(r[cols],col=="red")
fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
prefix <- c("Sentinel variant","CHR","POS","Mapped gene","Target","-log10(p)")
postfix <- c("</br>")
d1 <- list(x=d3c$x, y=d3c$y, z=d3c$log10p, text=as.list(apply(sapply(1:6,fixes,d3c),1,paste,collapse=" ")),
           type="scatter3d", mode="markers", name="cis")
d1$marker$symbol <- head(xdata[[1]]$marker$symbol,nrow(d3c))
d2 <- list(x=d3t$x, y=d3t$y, z=d3t$log10p, text=as.list(apply(sapply(1:6,fixes,d3t),1,paste,collapse=" ")),
           type="scatter3d", mode="markers", name="trans")
d2$marker$symbol <- head(xdata[[2]]$marker$symbol,nrow(d3t))
x <- SomaLogic$x
x$data[[1]] <- d1
x$data[[2]] <- d2
Olink <- SomaLogic
Olink$x <- x
sink("js.json")
toJSON(Olink,auto_unbox=TRUE,pretty=TRUE)
sink()

# sed -i 's|<\\/br>|\\u003c/br>|g' js.json


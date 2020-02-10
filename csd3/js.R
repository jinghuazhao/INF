# 10-2-2020 JHZ

library(jsonlite)
SomaLogic <- read_json("js.SomaLogic")
library(gap)
d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
r <-  mhtplot2d(d)
r <- within(r,{x=x/1.3e8; y=y/1.3e8; log10p <- abs(log10p)})
prefix <- c("Sentinel variant","CHR","POS","Mapped gene","Target","-log10(p)")
postfix <- c("</br>")
fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
cols <- c("id","chr1","x","gene","target","log10p","chr2","y","col")
Olink <- SomaLogic
a <- subset(r[cols],col=="red")
Olink$x$data[[2]] <- list(x=a$x, y=a$y, z=a$log10p, text=as.list(apply(sapply(1:6,fixes,a),1,paste,collapse=" ")),
                     type="scatter3d", mode="markers", name="cis")
Olink$x$data[[2]]$marker$symbol <- rep('circle',nrow(a))
Olink$x$data[[2]]$marker$size <- 3
b <- subset(r[cols],col=="blue")
Olink$x$data[[1]] <- list(x=b$x, y=b$y, z=b$log10p, text=as.list(apply(sapply(1:6,fixes,b),1,paste,collapse=" ")),
                     type="scatter3d", mode="markers", name="trans")
Olink$x$data[[1]]$marker$symbol <- rep('circle',nrow(b))
Olink$x$data[[1]]$marker$size <- 3
sink("js.json")
toJSON(Olink,auto_unbox=TRUE,pretty=TRUE)
sink()

# sed -i 's|<\\/br>|\\u003c/br>|g' js.json

# 10-2-2020 JHZ

library(jsonlite)
SomaLogic <- read_json("js.SomaLogic")
library(gap)
d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
r <-  mhtplot2d(d)
r <- within(r,{x=x/1.3e8; y=y/1.3e8; log10p <- abs(log10p)})
x <- with(SomaLogic,x)
prefix <- c("Sentinel variant","CHR","POS","Mapped gene","Target","-log10(p)")
postfix <- c("</br>")
fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
cols <- c("id","chr1","x","gene","target","log10p","chr2","y","col")
a <- subset(r[cols],col=="blue")
b <- subset(r[cols],col=="red")
x$data[[1]] <- list(x=a$x, y=a$y, z=a$log10p, text=as.list(apply(sapply(1:6,fixes,a),1,paste,collapse=" ")),
           type="scatter3d", mode="markers", name="trans")
x$data[[1]]$marker$symbol <- rep('circle',nrow(a))
x$data[[1]]$marker$size <- 3
x$data[[2]] <- list(x=b$x, y=b$y, z=b$log10p, text=as.list(apply(sapply(1:6,fixes,b),1,paste,collapse=" ")),
           type="scatter3d", mode="markers", name="cis")
x$data[[2]]$marker$symbol <- rep('circle',nrow(b))
x$data[[2]]$marker$size <- 3
Olink <- SomaLogic
Olink$x <- x
sink("js.json")
toJSON(Olink,auto_unbox=TRUE,pretty=TRUE)
sink()

# sed -i 's|<\\/br>|\\u003c/br>|g' js.json

# 9-2-2020 JHZ

R --no-save <<END
  library(gap)
  d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  r <-  mhtplot2d(d)
  r <- within(r,{log10p <- abs(log10p)})
  library(jsonlite)
  SomaLogic <- read_json("js.SomaLogic")
  xdata <- SomaLogic$x$data
  cols <- c("id","chr1","pos1","gene","target","log10p","chr2","pos2","col")
  d3c <- subset(r[cols],col=="blue")
  d3t <- subset(r[cols],col=="red")
  fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
  prefix <- c("Sentinel variant","CHR","POS","Mapped gene","Target","-log10(p)")
  postfix <- c("</br>")
  d1 <- list(x=d3c$pos1, y=d3c$pos2, z=d3c$log10p, text=as.list(apply(sapply(1:6,fixes,d3c),1,paste,collapse=" ")),
             type="scatter3d", mode="markers", name="cis")
  d1$marker$symbol <- head(xdata[[1]]$marker$symbol,nrow(d3c))
  d2 <- list(x=d3t$pos1, y=d3t$pos2, z=d3t$log10p, text=as.list(apply(sapply(1:6,fixes,d3t),1,paste,collapse=" ")),
             type="scatter3d", mode="markers", name="trans")
  d2$marker$symbol <- head(xdata[[2]]$marker$symbol,nrow(d3t))
  x <- SomaLogic$x
  x$data[[1]] <- d1
  x$data[[2]] <- d2
  Olink <- SomaLogic
  Olink$x <- x
  sink("js.json")
  toJSON(Olink,auto_unbox=TRUE)
  sink()
END
sed -i 's|<\\/br>|\\u003c/br>|g' js.json

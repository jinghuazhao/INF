# 7-2-2020 JHZ

R --no-save <<END
  library(gap)
  d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  r <-  mhtplot2d(d)
  r <- within(r,{log10p <- abs(log10p)})
  write.table(r,file="INF1.merge.2d",quote=FALSE,row.names=FALSE,sep=",")
END
sed -i 's/red/1/g;s/blue/2/g' INF1.merge.2d
R --no-save <<END
  library(jsonlite)
  SomaLogic <- read_json("js.SomaLogic")
  x <- SomaLogic$x
  data <- x$data
  data1 <- data[[1]]
  data2 <- data[[2]]
  Olink <- SomaLogic
  d3 <- read.csv("INF1.merge.2d",as.is=TRUE)
  prefix <- c("Sentinel variant","CHR","POS","Mapped gene","Target","-log10(p)")
  cols <- c("id","chr1","pos1","gene","target","log10p","chr2","pos2","col")
  postfix <- c("</br>")
  d3cis <- subset(d3[cols],col==1)
  d3trans <- subset(d3[cols],col==2)
  fixes <- function(col,d) paste(paste(prefix[col],d[,col],sep=":"),postfix)
  text1 <- as.list(apply(sapply(1:6,fixes,d3cis),1,paste,collapse=" "))
  text2 <- as.list(apply(sapply(1:6,fixes,d3trans),1,paste,collapse=" "))
  symbol1 <- head(data1$marker$symbol,nrow(d3cis))
  symbol2 <- head(data2$marker$symbol,nrow(d3trans))
  x1 <- with(d3cis, list(x=pos1, y=pos2, z=log10p))
  x2 <- with(d3trans, list(x=pos1, y=pos2, z=log10p))
  d1 <- within(data1,{x=x1; name="cis"; marker$symbol=symbol1; text=text1})
  d2 <- within(data2,{x <- x2; name <- "trans"; marker$symbol <- symbol2; text=text2})
  Olink$x$data <- list(d1,d2)
  sink("js.txt")
  toJSON(Olink)
  sink()
END

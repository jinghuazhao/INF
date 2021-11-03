LTBR <- function(xy,out="p2d.png")
{
  data <- data.frame(xy)
  x <- data[[1]]
  png(out,res=300,width=18,height=10,units="in")
  par(mfrow=c(2,1))
  plot(data[[2]]~x,ylim=c(-0.5,0.5),col="blue",xlab="Position",ylab="effect size",pch=15,axes=FALSE)
  points(data[[3]]~x,col="orange",pch=17)
  points(data[[4]]~x,col="red",pch=19)
  title("MS-LTBR-TNFB")
  legend("topleft",
         legend = c("MS","LTBR","TNFB"),
         col = c("blue","orange","red"),
         pch = c(15,17,19),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = FALSE,
         inset = c(0.1, 0.1))
  axis(1,at=x,tick=TRUE)
  axis(2)
  plot(data[[2]]~x,ylim=c(-0.5,0.5),col="blue",type="b",xlab="Position",ylab="effect size",pch=15,axes=FALSE)
  lines(data[[3]]~x,col="orange",pch=17)
  lines(data[[4]]~x,col="red",pch=19)
  title("MS-LTBR-TNFB")
  legend("topleft",
         legend = c("MS","LTBR","TNFB"),
         col = c("blue","orange","red"),
         pch = c(15,17,19),
         bty = "n",
         pt.cex = 2,
         cex = 1.2,
         text.col = "black",
         horiz = FALSE,
         inset = c(0.1, 0.1))
  axis(1,at=x,tick=TRUE)
  axis(2)
  dev.off()
}

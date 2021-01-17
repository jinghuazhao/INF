# 17-1-2021 JHZ

zld <- function(z)
{
  id <- with(z,index)
  rank <- 1:with(z,length(index))
  ldt <- ld[id,id]
  ldt[upper.tri(ldt, diag=TRUE)] <- NA
  cbind(rank,z[rank,c("index","rsid","z","log10bf","group","corr_group","prob_group","log10bf_group")],ldt)
}

pip.log10bf.plot <- function()
{
  png(paste0(pr,".png"),height=12,width=8,units="in",res=300)
  par(mfrow=c(2,1))
  with(snp, {
    plot(position/1e6,prob,cex=0.3,xlab="",ylab="PIP",axes=FALSE)
    points(position[prob>0.8]/1e6,prob[prob>0.8],cex=0.3,col="red")
    axis(2)
    plot(position/1e6,log10bf,cex=0.3,xlab="Position (MB)", ylab="log10(BF)", axes=FALSE)
    axis(1)
    axis(2)
  })
  dev.off()
}

options(digits=3, scipen=20, width=500)
pr <- Sys.getenv("pr")
z <- read.table(paste0(pr, ".z"), as.is=TRUE, header=TRUE)
ld <- read.table(paste0(pr, ".ld.gz"),col.names=with(z,rsid))
snp <- read.table(paste0(pr, ".snp"), as.is=TRUE, header=TRUE)
config <- read.table(paste0(pr,".config"),as.is=TRUE,header=TRUE,nrow=50)
if (file.exists(paste0(pr,".cred"))) cred <- read.table(paste0(pr,".cred"),as.is=TRUE,header=TRUE)

library(openxlsx)
xlsx <- paste0(pr, "-finemap.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
# snp
  d <- within(snp,{log10p_incl <- gap::log10p(mean_incl/sd_incl)})
  name_snp <- d[,setdiff(names(d),c("chromosome","position","allele1","allele2","maf","beta","se","rank"))]
  addWorksheet(wb, "snp", zoom=150)
  ord <- with(name_snp,order(-prob_group))
  writeDataTable(wb, "snp", name_snp[ord,])
# pip.plot
  pip.log10bf.plot()
  addWorksheet(wb, "pip.plot", zoom=150)
  insertImage(wb, "pip.plot", paste0(pr,".png"),height=12,width=8)
# config
  addWorksheet(wb, "config", zoom=150)
  writeDataTable(wb, "config", config)
# cred
  if (exists("cred")) {
    addWorksheet(wb, "cred", zoom=150)
    writeDataTable(wb, "cred", cred)
  }
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

topz <- subset(snp, abs(z)>=3.290527)
if (nrow(topz) > 0)
{
# topz
  zldz <- zld(topz)
  topsnp <- head(snp, nrow(topz))
  zldsnp <- zld(topsnp)
  xlsx <- paste0(pr, "-zld.xlsx")
  unlink(xlsx, recursive = FALSE, force = TRUE)
  wb <- createWorkbook(xlsx)
  addWorksheet(wb, "topz", zoom=150)
  writeDataTable(wb, "topz", zldz)
# topsnp
  addWorksheet(wb, "topsnp", zoom=150)
  writeDataTable(wb, "topsnp", zldsnp)
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}

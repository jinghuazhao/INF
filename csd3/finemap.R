# 17-9-2019 JHZ

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
topz <- subset(snp, abs(z)>=6.47)
zldz <- zld(topz)
topsnp <- head(snp, nrow(topz))
zldsnp <- zld(topsnp)
pip.log10bf.plot()
load(paste0(pr,".rda"))
snp <- within(snp, {rank <- 1:nrow(snp)})
d <- merge(snpid_rsid,snp,by="rsid",all.y=TRUE)
ord <-order(with(d,rank))
snp <- d[ord,setdiff(names(d),c("chromosome","position","allele1","allele2","maf","beta","se","rank"))]
config <- read.table(paste0(pr,".config"),as.is=TRUE,header=TRUE,nrows=31)
if (file.exists(paste0(pr,".cred"))) cred <- read.table(paste0(pr,".cred"),as.is=TRUE,header=TRUE)

library(openxlsx)
xlsx <- paste0(pr, "-finemap.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, "snp")
writeDataTable(wb, "snp", snp)
addWorksheet(wb, "pip.plot")
insertImage(wb, "pip.plot", paste0(pr,".png"),height=12,width=8)
addWorksheet(wb, "config")
writeDataTable(wb, "config", config)
if (exists("cred")) {
  addWorksheet(wb, "cred")
  writeDataTable(wb, "cred", cred)
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

xlsx <- paste0(pr, ".zld.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, "topz")
writeDataTable(wb, "topz", zldz)
addWorksheet(wb, "topsnp")
writeDataTable(wb, "topsnp", zldsnp)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

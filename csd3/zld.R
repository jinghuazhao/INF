# 15-9-2019 JHZ

zld <- function(z)
{
  id <- with(z,index)
  rank <- 1:with(z,length(index))
  ldt <- ld[id,id]
  ldt[upper.tri(ldt, diag=TRUE)] <- NA
  cbind(rank,z[rank,c("index","rsid","z","log10bf","group","corr_group","prob_group","log10bf_group")],ldt)
}
options(digits=3, scipen=20, width=10000)
pr <- Sys.getenv("pr")
z <- read.table(paste0(pr, ".z"), as.is=TRUE, header=TRUE)
ld <- read.table(paste0(pr, ".ld"),col.names=with(z,rsid))
snp <- read.table(paste0(pr, ".snp"), as.is=TRUE, header=TRUE)
topz <- subset(snp, abs(z)>=6.47)
zldz <- zld(topz)
topsnp <- head(snp, nrow(topz))
zldsnp <- zld(topsnp)
library(openxlsx)
xlsx <- paste0(pr, ".zld.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
f <- make.names(pr)
addWorksheet(wb, paste0(f, ".topz"))
writeDataTable(wb, paste0(f, ".topz"), zldz)
addWorksheet(wb, paste0(f, ".topsnp"))
writeDataTable(wb, paste0(f, ".topsnp"), zldsnp)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

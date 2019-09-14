# 14-9-2019 JHZ

options(digits=3, scipen=20, width=10000)
pr <- Sys.getenv("pr")
z <- read.table(paste0(pr, ".z", as.is=TRUE, header=TRUE)
ld <- read.table(paste0(pr, ".ld"),col.names=with(z,rsid))
snp <- read.table(paste0(pr, ".snp"), as.is=TRUE, header=TRUE)
z <- subset(snp, abs(z) > 6.47)
id <- with(z,index)
rank <- 1:with(z,length(index))
ldt <- ld[id,id][rank,rank]
ldt[upper.tri(ldt, diag=TRUE)] <- NA
chk <- cbind(rank,z[rank,c("index","rsid","z","group","corr_group","prob_group")],ldt)
library(openxlsx)
xlsx <- paste0(pr, "-zld.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
f <- make.names(pr)
addWorksheet(wb, paste0(f, ".zld"))
writeDataTable(wb, paste0(f, ".zld"), chk)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

# 10-9-2019 JHZ

options(digits=3, scipen=20, width=10000)
pr <- Sys.getenv("pr")
cat(pr, "\n")
ld <- read.table(paste0(pr, ".ld"))
snp <- read.table(paste0(pr, ".snp"), as.is=TRUE, header=TRUE)
z <- subset(snp, abs(z) > 6.47)
id <- as.integer(row.names(z))
ldt <- ld[id,id]
ldt[upper.tri(ldt, diag=TRUE)] <- NA
colnames(ldt) <- with(z, rsid)
chk <- cbind(z[c("rsid","z")], ldt)
library(openxlsx)
xlsx <- paste0(pr, "-zld.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
f <- make.names(pr)
addWorksheet(wb, paste0(f, ".zld"))
writeDataTable(wb, paste0(f, ".zld"), chk)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)


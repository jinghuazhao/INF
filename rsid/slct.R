# 6-5-2020 JHZ

options(scipen=20, width=2000)
pr <- Sys.getenv("pr")
jma <- read.delim(paste0(pr,".jma.cojo"),as.is=TRUE)
ldr <- read.delim(paste0(pr,".ldr.cojo"), as.is=TRUE)
tbl <- jma[setdiff(names(jma),c("b","se","p"))]
cred <- gap::cs(tbl, b="bJ", se="bJ_se", cutoff=0.95)
require(openxlsx)
xlsx <- paste0(pr,"-slct.xlsx")
wb <- createWorkbook(xlsx)
f <- make.names(pr)
addWorksheet(wb, "jma")
writeDataTable(wb, "jma", jma)
addWorksheet(wb, "ldr")
writeDataTable(wb, "ldr", ldr)
addWorksheet(wb, "cs")
writeDataTable(wb, "cs",  cred)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

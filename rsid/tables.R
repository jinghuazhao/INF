options(scipen=20, width=2000)
require(openxlsx)

url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
summary <-
studies <- 
inf1 <-
pqtls <-
cojo <-
eqtls <-
eqtlstudies <-
transpqtls <-
knownpqtls <-
pqtlstudies <-
pqtldisease <-
efo <- subset(read.xlsx(url, sheet="EFO", colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:4), rows=c(1:78)),!is.na(MRBASEID))

prefix <- "ST"
titles <- c("summary","study information","panel information",
            "pQTL summary","conditional analysis",
            "eQTL overlap","eQTL studies","trans pQTL mapping","transpqtls","previous pQTLs","previous pQTL studies","SomaLogic replication",
            "Disease GWAS overlap")
sheetnames <- make.names(paste(prefix,titles,sep=" - "))

xlsx <- paste0("scallop.xlsx")
wb <- createWorkbook(xlsx)
for (sheet in c("summary","studies","inf1","pqtls","cojo","eqtls","eqtlstudies","transpqtl","pqtlstudies","knownpqtls","pqtldisease"))
{
  addWorksheet(wb, sheet)
  writeDataTable(wb, sheetnames[1], get(sheet))
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

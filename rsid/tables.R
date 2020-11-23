options(scipen=20, width=2000)
require(openxlsx)

url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
sheets <- c("summary","studies","inf1","pqtls","cojo","eqtls","eqtlstudies","transpqtls","knownpqtls","pqtlstudies","somalogic","pqtldisease")
prefix <- "ST"
titles <- c("summary","study information","panel information",
            "pQTL summary","conditional analysis",
            "eQTL overlap","eQTL studies","trans pQTL mapping","previous pQTLs","previous pQTL studies","SomaLogic replication",
            "Disease GWAS overlap")
sheetnames <- paste(prefix,titles,sep=" - ")

summary <- data.frame(Sheet=sheets,Description=sheetnames)
studies <- read.xlsx(url, sheet="Studies", colNames=TRUE, skipEmptyRows=TRUE, cols=1:3, rows=2:15)
inf1 <- subset(read.xlsx(url, sheet="INF1", colNames=TRUE, skipEmptyRows=TRUE, cols=1:12, rows=1:93),uniprot!="P23560")
pqtls <- read.xlsx(url, sheet="pQTLs", colNames=TRUE, skipEmptyRows=TRUE, cols=1:21, rows=1:181)
cojo <- read.xlsx(url, sheet="cojo", colNames=TRUE, skipEmptyRows=TRUE, cols=1:22, rows=1:228)
eqtls <- subset(within(read.xlsx(url, sheet="eQTLs", colNames=TRUE, skipEmptyRows=TRUE, cols=1:36, rows=1:214),{allele1 <- a1;allele2 <- a2}),
                select=-c(a1,a2))
eqtlstudies <- data.frame(Studies="Studies",pmid="pmid")
transpqtls <- data.frame(Studies="Studies",pmid="pmid")
knownpqtls <- data.frame(Studies="Studies",pmid="pmid")
pqtlstudies <- data.frame(Studies="Studies",pmid="pmid")
pqtldisease <- data.frame(Studies="Studies",pmid="pmid")
somalogic <- data.frame(Studies="Studies",pmid="pmid")
efo <- subset(read.xlsx(url, sheet="EFO", colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:4), rows=c(1:78)),!is.na(MRBASEID))

xlsx <- paste0("work/SCALLOP-INF.xlsx")
wb <- createWorkbook(xlsx)
for (i in 1:length(sheets))
{
  cat(sheets[i],sheetnames[i],"\n")
  addWorksheet(wb, sheetnames[i])
  writeDataTable(wb, sheetnames[i], get(sheets[i]))
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

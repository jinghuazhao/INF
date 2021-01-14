#options(scipen=20, width=2000)

require(openxlsx)
url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
read.sheet <- function(url=url,sheet,cols,rows)
{
  read.xlsx(url,sheet=sheet,colNames=TRUE,cols=cols,rows=rows,skipEmptyRows=TRUE)
}

summary <- read.sheet(sheet="Summary", cols=1:2, rows=2:36)
inf1 <- subset(read.sheet(sheet="INF1", cols=1:12, rows=2:94),uniprot!="P23560")
studies <- read.sheet(sheet="Studies", cols=1:3, rows=2:14)
pqtls <- read.sheet(sheet="pQTLs", cols=1:21, rows=2:182)
interval <- read.sheet(sheet="INTERVAL", cols=1:12, rows=2:29)
otherknownloci <- read.sheet(sheet="OtherKnownLoci", cols=1:12, rows=2:102)
cvd1 <- read.sheet(sheet="CVD1", cols=1:12, rows=2:53)
aristotle <- read.sheet(sheet="ARISTOTLE", cols=1:14, rows=2:182)
cojo <- read.sheet(sheet="cojo", cols=1:19, rows=2:229)
h2pve <- read.sheet(sheet="h2pve", cols=1:10, rows=2:93)
vep <- read.sheet(sheet="VEP", cols=1:69, rows=2:261)
eqtls <- read.sheet(sheet="eQTLs", cols=1:24, rows=2:24)
reactome <- read.sheet(sheet="reactome", cols=1:19, rows=2:589)
great <- read.sheet(sheet="GREAT", cols=1:24, rows=2:30)
garfield <- read.sheet(sheet="GARFIELD", cols=1:18, rows=2:7037)
fusion <- read.sheet(sheet="FUSION", cols=1:26, rows=2:117)
efo <- subset(read.sheet(sheet="EFO", cols=c(1:4), rows=c(2:79)),!is.na(MRBASEID))
smr <- read.sheet(sheet="SMR", cols=1:27, rows=2:83)
pqtlmr_immue <- read.sheet(sheet="pqtlMR-immune", cols=1:7, rows=2:67)
pqtlmr_misc <- read.sheet(sheet="pqtlMR-misc", cols=1:7, rows=2:39)
ivw <- read.sheet(sheet="IVW", cols=1:7, rows=2:12)
gsmr <- read.sheet(sheet="GSMR", cols=1:12, rows=2:55)
crp <- read.sheet(sheet="CRP", cols=1:15, rows=2:30)
mrc2 <- read.sheet(sheet="MRC2", cols=1:7, rows=2:8)
mvmr <- read.sheet(sheet="MVMR", cols=1:9, rows=2:6)
genedrugbank <- read.sheet(sheet="geneDrugbank", cols=1:7, rows=2:72)
eqtlstudies <- data.frame(Studies="Studies",pmid="pmid")
transpqtls <- data.frame(Studies="Studies",pmid="pmid")
knownpqtls <- data.frame(Studies="Studies",pmid="pmid")
pqtlstudies <- data.frame(Studies="Studies",pmid="pmid")
pqtldisease <- data.frame(Studies="Studies",pmid="pmid")
somalogic <- data.frame(Studies="Studies",pmid="pmid")

xlsx <- paste0("work/SCALLOP-INF.xlsx")
outsheets <- c("summary","studies","inf1","pqtls","cojo","eqtls","eqtlstudies","knownpqtls","pqtlstudies","somalogic","pqtldisease")
prefix <- "ST"
titles <- c("summary","study information","panel information",
            "pQTL summary","conditional analysis",
            "eQTL overlap","eQTL studies","trans pQTL mapping","previous pQTLs","previous pQTL studies","SomaLogic replication",
            "Disease GWAS overlap")
sheetnames <- paste(prefix,titles,sep=" - ")
wb <- createWorkbook(xlsx)
for (i in 1:length(sheets))
{
  cat(sheets[i],sheetnames[i],"\n")
  addWorksheet(wb, sheetnames[i])
  writeDataTable(wb, sheetnames[i], get(sheets[i]))
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

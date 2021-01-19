#options(scipen=20, width=2000)
require(openxlsx)
url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
read.sheet <- function(sheet,cols,rows) read.xlsx(url,sheet=sheet,colNames=TRUE,cols=cols,rows=rows,skipEmptyRows=TRUE)

 summary <- read.sheet("Summary", 1:2, 2:36)
    inf1 <- subset(read.sheet("INF1", 1:12, 2:94),uniprot!="P23560")
 studies <- read.sheet("Studies", 1:3, 2:14)
   pqtls <- read.sheet("pQTLs", 1:21, 2:182)
interval <- read.sheet("INTERVAL", 1:12, 2:29)
      os <- read.sheet("OtherStudies", 1:12, 2:102)
    cvd1 <- read.sheet("CVD1", 1:12, 2:53)
aristotl <- read.sheet("ARISTOTLE", 1:14, 2:182)
    cojo <- read.sheet("cojo", 1:19, 2:229)
   h2pve <- read.sheet("h2pve", 1:10, 2:93)
     vep <- read.sheet("VEP", 1:69, 2:261)
   eqtls <- read.sheet("eQTLs", 1:24, 2:24)
reactome <- read.sheet("Reactome", 1:19, 2:589)
   great <- read.sheet("GREAT", 1:24, 2:30)
garfield <- read.sheet("GARFIELD", 1:18, 2:7037)
  fusion <- read.sheet("FUSION", 1:26, 2:117)
     efo <- subset(read.sheet("EFO", 1:4, 2:79),!is.na(MRBASEID))
     smr <- read.sheet("SMR", 1:27, 2:83)
mr_immue <- read.sheet("pqtlMR-immune", 1:7, 2:67)
 mr_misc <- read.sheet("pqtlMR-misc", 1:7, 2:39)
     ivw <- read.sheet("IVW", 1:7, 2:12)
    gsmr <- read.sheet("GSMR", 1:12, 2:55)
     crp <- read.sheet("CRP", 1:15, 2:30)
    mrc2 <- read.sheet("MRC2", 1:7, 2:8)
    mvmr <- read.sheet("MVMR", 1:9, 2:6)
     gdb <- read.sheet("geneDrugbank", 1:7, 2:72)

knownpqtls_dup <- within(rbind(interval,os,cvd1),{
                         Sentinels <- gsub(" ","",Sentinels)
                         UniProt <- gsub(" ","",UniProt)
                         Protein <- gsub(" ","",Protein)
                         SNPid <- gsub(" ","",SNPid)
                         Source <- gsub(" ","",Source)
                  })
knownpqtls <- unique(knownpqtls_dup[c("Sentinels","SNPid","UniProt","Protein")])
pqtlstudies <- unique(knownpqtls_dup[c("Source","PMID")])
ord <- with(pqtlstudies,order(PMID))
pqtlstudies <- pqtlstudies[ord,]
rownames(pqtlstudies) <- seq(nrow(pqtlstudies))

url <- "https://jhz22.user.srcf.net/pqtl-immune_infection_edited.xlsx"
  metal <- read.sheet("METAL",1:20,1:181)
  short <- read.sheet("short",1:51,1:220)

pqtldisease <- subset(short,Keep==1,select=c(MarkerName,nprots,prots,Allele1,Allele2,Effects,SEs,cistrans,trait,efo,study,pmid,dataset,infection))

outsheets <- c("summary","studies","inf1","pqtls","cojo","knownpqtls","pqtlstudies","interval","eqtls","pqtldisease")
titles <- c("summary","study information","panel information","pQTLs","conditional analysis",
            "known pQTLs","previous pQTL studies","SomaLogic replication","eQTLs","Disease GWAS overlap")
description=paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))
uppered <- c("PQTLs","EQTLs")
description[description%in%uppered] <- titles[description%in%uppered]
summary <- data.frame(Sheetnames=paste0("ST",1:length(outsheets)),Description=description)
xlsx <- paste0("work/SCALLOP-INF.xlsx")
wb <- createWorkbook(xlsx)
options("openxlsx.borderColour"="#4F80BD")
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
for (i in 1:length(outsheets))
{
  sheetnames <- with(summary,paste0(Sheetnames,"-",Description)[i])
  cat(outsheets[i],sheetnames,"\n")
  addWorksheet(wb, sheetnames, zoom=150)
  writeData(wb, sheetnames, sheetnames, startCol=1, startRow=1)
  writeDataTable(wb, sheetnames, get(outsheets[i]), startCol=1, startRow=2,
                 headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
}
sheets(wb)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

novelpqtls <- subset(pqtls,!paste0(prot,"-",rsid)%in%with(knownpqtls,paste0(Protein,"-",Sentinels)),select=c(MarkerName,rsid,prot,uniprot))
ord <- with(novelpqtls,order(prot,MarkerName))
write.xlsx(cbind(no=1:nrow(novelpqtls),novelpqtls[ord,]), file="novelpqtls.xlsx", colNames=TRUE,
           borders="surrounding", headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")

options(width=2000)
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
mr_immun <- read.sheet("pqtlMR-immune", 1:7, 2:67)
 mr_misc <- read.sheet("pqtlMR-misc", 1:7, 2:39)
     ivw <- read.sheet("IVW", 1:8, 2:19)
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

INF <- Sys.getenv("INF")
options("openxlsx.borderColour"="#4F80BD")
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
novelpqtls <- subset(pqtls,!paste0(prot,"-",rsid)%in%with(knownpqtls,paste0(Protein,"-",Sentinels)),select=c(MarkerName,rsid,prot,uniprot))
ord <- with(novelpqtls,order(prot,MarkerName))
write.xlsx(cbind(no=1:nrow(novelpqtls),novelpqtls[ord,]), file=file.path(INF,"NG","novelpqtls.xlsx"), colNames=TRUE,
           borders="surrounding", headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")

url <- "https://jhz22.user.srcf.net/pqtl-immune_infection_edited.xlsx"
metal <- read.sheet("METAL",1:20,1:181)
short <- read.sheet("short",1:51,1:220)
pqtldisease <- subset(short,Keep==1,select=c(MarkerName,nprots,prots,Allele1,Allele2,Effects,SEs,cistrans,trait,efo,study,pmid,dataset,infection))
credibleset <- read.table(file.path(INF,"work","INF1.merge-rsid.cs"),col.names=c("prot","MarkerName","CredibleSet"),sep="\t")
pqtls <- merge(pqtls,credibleset,by.x=c("prot","rsid"),by.y=c("prot","MarkerName"))
coloc <- read.delim(file.path(INF,"coloc","GTEx.tsv"))
cs95 <- read.delim(file.path(INF,"coloc","cis-eQTL_table.tsv"))
cs95 <- data.frame(rsidProt=rownames(cs95),cs95)

outsheets <- c("summary","studies","inf1","interval","os","cvd1","aristotl",
               "pqtls","cojo","knownpqtls","pqtlstudies","coloc","cs95","smr","pqtldisease",
               "vep","garfield",
               "mr_immun","ivw","mr_misc","gsmr","gdb")
titles <- c("summary","study information","panel information","INTERVAL study","Other studies","SCALLOP-CVD1","ARISTOTLE study",
            "pQTLs","conditional analysis",
            "known pQTLs","previous pQTL studies","GTEx coloc","GTEx coloc 95%CS","SMR","Disease GWAS overlap",
            "VEP annotation","GARFIELD outputs",
            "pQTL-immune-MR","IVW","pQTL-misc-MR","GSMR-FEV1CVD","geneDrugbank")
description=paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))
uppered <- c("PQTLs")
description[description%in%uppered] <- titles[description%in%uppered]
n0 <- 7
prefix <- c(paste0(toupper(substr(outsheets, 1, 1)), substr(outsheets, 2, nchar(outsheets)))[1:n0],paste0("ST",1:(length(outsheets)-n0)))
summary <- data.frame(Sheetnames=prefix,Description=description)
xlsx <- file.path(INF,"NG","SCALLOP-INF.xlsx")
wb <- createWorkbook(xlsx)
for (i in 1:length(outsheets))
{
  sheetnames <- with(summary[i,], ifelse(i<=n0, Description, paste0(Sheetnames,"-",Description)))
  cat(sheetnames,"\n")
  addWorksheet(wb, sheetnames, zoom=150)
  writeData(wb, sheetnames, sheetnames, startCol=1, startRow=1)
  writeDataTable(wb, sheetnames, get(outsheets[i]), startCol=1, startRow=2,
                 headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
}
data.frame(sheets(wb))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

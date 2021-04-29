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
     vep <- read.sheet("VEP", 1:27, 2:182)
   eqtls <- read.sheet("eQTLs", 1:24, 2:24)
reactome <- read.sheet("Reactome", 1:19, 2:589)
   great <- read.sheet("GREAT", 1:24, 2:101)
  great3 <- read.sheet("IL12B-KITLG-TNFSF10", 1:25, 2:38)
garfield <- read.sheet("GARFIELD", 1:18, 2:7037); o <- with(garfield, order(Pvalue)); garfield <- garfield[o,]
  fusion <- read.sheet("FUSION", 1:26, 2:117)
     efo <- subset(read.sheet("EFO", 1:4, 2:79),!is.na(MRBASEID))
     smr <- read.sheet("SMR", 1:27, 2:83)
    gsmr <- read.sheet("GSMR", 1:12, 2:55)
     crp <- read.sheet("CRP", 1:15, 2:30)
     gdb <- read.sheet("geneDrugbank", 1:7, 2:72)
     at1 <- readWorkbook(xlsxFile=url,sheet="Annotrans1"); #names(at1) <- replace(names(at1),grepl("^[X]",names(at1)),"")
     at2 <- readWorkbook(xlsxFile=url,sheet="Annotrans2"); #names(at2) <- replace(names(at2),grepl("^[X]",names(at2)),"")
     at3 <- readWorkbook(xlsxFile=url,sheet="Annotrans3"); #names(at3) <- replace(names(at3),grepl("^[X]",names(at3)),"")

INF <- Sys.getenv("INF")
read_table <- function(f, exprs="pval <= 0.05/nrow(t)")
{
   t <- within(read.delim(f), {Flag=""})
   o <- with(t,order(outcome,pval))
   cond <- str2expression(exprs)
   n <- subset(t, eval(cond))
   x <- with(t,eval(cond))
   t[x, "Flag"] <- "x"
   t[o,]
}
mr_immun <- read_table(file.path(INF,"mr","pQTLs","pQTL-efo.txt"))
mr_misc <- read_table(file.path(INF,"mr","pQTLs","pQTL-ieu-FEV1.txt"))
mr <- read_table(file.path(INF,"mr","efo-result.txt"),exprs="cistrans!=\"pan\" & pval <= 0.05/nrow(t)")

pav <- merge(within(pqtls,{prot_rsid=paste0(prot,"-",rsid)}),within(vep,{prot_rsid=paste0(Protein,"-",vep[["#Uploaded_variation"]])}),by="prot_rsid") 
data.frame(table(subset(pav,cis.trans=="cis")$Consequence))

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

options("openxlsx.borderColour"="#4F80BD")
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
novelpqtls <- subset(within(pqtls,{
                                   chrpos=paste0(Chromosome,":",Position)
                                   a1a2=paste0(toupper(Allele1),"/",toupper(Allele2))
                                   bse=paste0(round(Effect,3)," (",round(StdErr,3),")")
                                   log10p=-log.P.
                                  }),
                     !paste0(prot,"-",rsid)%in%with(knownpqtls,paste0(Protein,"-",Sentinels)),
                     select=c(prot,uniprot,chrpos,rsid,a1a2,bse,log10p))
ord <- with(novelpqtls,order(prot,chrpos))
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
               "pqtls","cojo","knownpqtls","pqtlstudies","smr","coloc","cs95","pqtldisease",
               "vep","great3","garfield",
               "mr_immun","mr","mr_misc","gsmr","gdb","at1","at2","at3","reactome","great","efo")
titles <- c("summary","study information","panel information","INTERVAL study","Other studies","SCALLOP-CVD1","ARISTOTLE study",
            "pQTLs","conditional analysis",
            "known pQTLs","previous pQTL studies","SMR","GTEx coloc","GTEx coloc 95%CS","Disease GWAS overlap",
            "VEP annotation","IL12B-KITLG-TNFSF10","GARFIELD outputs",
            "pQTL-immune-MR","MR results","pQTL-misc-MR","GSMR-FEV1CVD","geneDrugbank","Annotrans-1","Annotrans-2","Annotrans-3","Reactome","GREAT","EFO")
description=paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))
uppered <- c("PQTLs")
description[description%in%uppered] <- titles[description%in%uppered]
n0 <- 7
n1 <- 16
prefix <- c(paste0(toupper(substr(outsheets, 1, 1)), substr(outsheets, 2, nchar(outsheets)))[1:n0],
            paste0("ST",1:n1),
            paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))[(n0+n1+1):length(outsheets)]
          )
summary <- data.frame(Sheetnames=prefix,Description=description)
xlsx <- file.path(INF,"NG","SCALLOP-INF.xlsx")
wb <- createWorkbook(xlsx)
for (i in 1:length(outsheets))
{
  sheetnames <- with(summary[i,], ifelse(i<=n0|i>n0+n1, Description, paste0(Sheetnames,"-",Description)))
  cat(sheetnames,"\n")
  if (i<=n0+n1 | i>n0+n1+3)
  {
    addWorksheet(wb, sheetnames, zoom=150)
    writeData(wb, sheetnames, sheetnames, xy=c(1,1),
                  headerStyle=createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
    body <- get(outsheets[i])
    writeDataTable(wb, sheetnames, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
    freezePane(wb, sheetnames, firstCol=TRUE, firstActiveRow=3)
    width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
  # width_vec_header <- nchar(colnames(body))+2
    setColWidths(wb, sheetnames, cols = 1:ncol(body), widths = width_vec)
    writeData(wb, sheetnames, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
  } else {
    sheet <- paste0("Annotrans-",i-n0-n1)
    addWorksheet(wb,sheet,gridLines=FALSE)
    writeData(wb,sheet,paste0("Trans-pQTL annotation-",i-n0-n1), xy=c(1,1),
                 headerStyle=createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
    writeDataTable(wb,sheet, get(paste0("at",i-n0-n1)), xy=c(1,2), firstColumn=TRUE, bandedRows=FALSE)
  }
}
sheets_wb <- sheets(wb)
data.frame(sheets_wb)

bStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
hStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
conditionalFormatting(wb, sheets_wb[grepl("GARFIELD",sheets_wb)], cols = 4, rows = 3:nrow(garfield), rule = "<=1e-5", style = hStyle)
conditionalFormatting(wb, sheets_wb[grepl("immune-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_immun), rule = "==\"x\"", style = hStyle)
conditionalFormatting(wb, sheets_wb[grepl("MR results",sheets_wb)], cols = 9, rows = 3:nrow(mr), rule = "==\"x\"", style = hStyle)
conditionalFormatting(wb, sheets_wb[grepl("misc-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_misc), rule = "==\"x\"", style = hStyle)

saveWorkbook(wb, file=xlsx, overwrite=TRUE)

# mr_immun <- read.sheet("pqtlMR-immune", 1:7, 2:67)
#  mr_misc <- read.sheet("pqtlMR-misc", 1:7, 2:39)
#      ivw <- read.sheet("IVW", 1:8, 2:19)
#     mrc2 <- read.sheet("MRC2", 1:7, 2:8)
#     mvmr <- read.sheet("MVMR", 1:9, 2:6)

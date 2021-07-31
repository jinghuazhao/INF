options(width=2000)
require(openxlsx)
url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
INF <- Sys.getenv("INF")
url <- file.path(INF,"work","INF1.latest.xlsx")
read.sheet <- function(sheet,cols,rows) read.xlsx(url,sheet=sheet,colNames=TRUE,cols=cols,rows=rows,skipEmptyRows=TRUE)
require(dplyr)
require(stringr)
gap_inf1 <- gap.datasets::inf1[c("uniprot", "prot", "target.short")]
 summary <- read.sheet("Summary", 1:2, 2:36)
    inf1 <- subset(read.sheet("INF1", 1:12, 2:94),uniprot!="P23560") %>% select(-panel); names(inf1)[2] <- "Protein"
 studies <- read.sheet("Studies", 1:3, 2:15)
   pqtls <- merge(read.sheet("pQTLs", 1:21, 2:182),gap_inf1[c("prot","target.short")],by="prot")
interval <- merge(within(read.sheet("INTERVAL", 1:12, 2:29),{Protein <- gsub(" ", "", Protein)}),
                  gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
      os <- merge(read.sheet("OtherStudies", 1:12, 2:102),gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
    cvd1 <- merge(read.sheet("CVD1", 1:12, 2:53),gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
aristotl <- merge(read.sheet("ARISTOTLE", 1:14, 2:182), gap_inf1[c("prot","target.short")], by.x="Protein", by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
    cojo <- merge(read.sheet("cojo", 1:19, 2:229),gap_inf1[c("prot","target.short")],by="prot") %>%
            mutate(prot=target.short) %>% rename(Protein=prot) %>% select(-target.short)
   h2pve <- read.sheet("h2pve", 1:10, 2:93)
     vep <- merge(read.sheet("VEP", 1:27, 2:182),gap_inf1,by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
   eqtls <- read.sheet("eQTLs", 1:24, 2:24)
reactome <- read.sheet("Reactome", 1:19, 2:589)
garfield <- read.sheet("GARFIELD", 1:18, 2:3017) %>%
            select(ID,PThresh,Pvalue,Annotation,Celltype,Tissue,Type,Category,OR,Beta,SE,CI95_lower,CI95_upper,NAnnotThesh,NAnnot,NThresh,N,linkID)
  fusion <- read.sheet("FUSION", 1:26, 2:117)
     smr <- merge(read.sheet("SMR", 1:27, 2:83),gap_inf1,by="prot") %>%
            mutate(prot=target.short) %>% rename(Protein=prot) %>% select(-target.short)
            d <- read.sheet("GSMR", 1:12, 2:55)
            na1 <- with(d,is.na(Exposure1))
            na2 <- with(d,is.na(Exposure2))
            d[na1,"Exposure1"] <- d[na1,"Exposure2"]
            d[na2,"Exposure2"] <- d[na2,"Exposure1"]
    gsmr <- merge(d, gap_inf1[c("prot","target.short")],by.x="Exposure1",by.y="prot") %>%
            mutate(Exposure1=target.short,Exposure2=target.short) %>% rename(Protein1=Exposure1,Protein2=Exposure2) %>%
            select(-target.short)
    gsmr_efo <- read.delim(file.path(INF,"mr","gsmr","out","5e-8","gsmr-efo.txt"))
    crp <- read.sheet("CRP", 1:15, 2:30)
    gdb <- read.sheet("geneDrugbank", 1:7, 2:72)
    at1 <- readWorkbook(xlsxFile=url,sheet="Annotrans1"); #names(at1) <- replace(names(at1),grepl("^[X]",names(at1)),"")
    at2 <- readWorkbook(xlsxFile=url,sheet="Annotrans2"); #names(at2) <- replace(names(at2),grepl("^[X]",names(at2)),"")
    at3 <- readWorkbook(xlsxFile=url,sheet="Annotrans3"); #names(at3) <- replace(names(at3),grepl("^[X]",names(at3)),"")

great3 <- read.delim(file.path(INF,"GREAT","IL12B-KITLG-TNFSF10.tsv")) %>%
          mutate(fdr=p.adjust(BinomP,method="fdr")) %>% arrange(fdr)
great <- read.delim(file.path(INF,"GREAT","cistrans.tsv")) %>%
          mutate(fdr=p.adjust(BinomP,method="fdr")) %>% arrange(fdr)

read_table <- function(f, exprs="pval <= 0.05/nrow(t)")
{
  addflag <- function(exprs)
  {
    t <- within(read.delim(f), {flag=""})
    x <- with(t,eval(str2expression(exprs)))
    x <- with(t,eval(str2expression(e)))
    t[x, "flag"] <- "x"
  }
  t <- within(read.delim(f), {fdr <- p.adjust(pval,method="fdr")})
# t <- addflag(exprs)
  t
}

mr_immun <- merge(read_table(file.path(INF,"mr","pQTLs","pQTL-efo.txt")),gap_inf1,by.x="exposure",by.y="prot") %>%
            mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)
mr_misc <- merge(read_table(file.path(INF,"mr","pQTLs","pQTL-ieu-FEV1.txt")),gap_inf1,by.x="exposure",by.y="prot") %>%
           mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)
mr <- merge(read_table(file.path(INF,"mr","efo-result.txt"),
            exprs="cistrans!=\"pan\" & pval <= 0.05/nrow(subset(t,cistrans!=\"pan\"))"),
            gap_inf1, by.x="exposure",by.y="prot") %>%
      mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)

protein_correlation <- read.csv(file.path(INF,"coffeeprot","table_complex.csv")) %>%
                       arrange(desc(cor)) %>%
                       mutate(varID1=toupper(varID1),varID2=toupper(varID2),VarVar=toupper(VarVar))
protein_dgi <- read.csv(file.path(INF,"coffeeprot","protein_annotated.csv")) %>%
               select(varID,ID,HPA_IF_protein_location,CP_loc,inDGIdb) %>%
               mutate(gene=toupper(ID)) %>%
               left_join(gap.datasets::inf1[c("gene","target.short")]) %>%
               rename(Protein=target.short) %>%
               select(-c(varID,ID,gene))
pqtl_annotation <- read.csv(file.path(INF,"coffeeprot","table_qtl_processed.csv")) %>%
                   rename(gene=gene_symbol) %>%
                   left_join(gap.datasets::inf1[c("gene","target.short")]) %>%
                   rename(Protein=target.short) %>%
                   select(-pvalue)
protein_dgi <- protein_dgi %>% select(Protein,names(protein_dgi))
pqtl_annotation <- pqtl_annotation %>% select(Protein,names(pqtl_annotation))
pqtl_impact <- read.csv(file.path(INF,"rsid","variant_effect_impact.csv"))

pav <- merge(within(pqtls,{prot_rsid=paste0(prot,"-",rsid)}),
             within(vep,{prot_rsid=paste0(Protein,"-",vep[["#Uploaded_variation"]])}),by="prot_rsid")
data.frame(table(subset(pav,cis.trans=="cis")$Consequence))

knownpqtls_dup <- rbind(interval,os,cvd1)
knownpqtls <- unique(knownpqtls_dup[c("Sentinels","SNPid","UniProt","Protein")])
pqtlstudies <- unique(knownpqtls_dup[c("Source","PMID")]) %>% arrange(PMID)
rownames(pqtlstudies) <- seq(nrow(pqtlstudies))

options("openxlsx.borderColour"="#4F80BD")
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
url <- "https://jhz22.user.srcf.net/pqtl-immune_infection_edited.xlsx"
url <- file.path(INF,"work","pqtl-immune_infection_edited.xlsx")
credibleset <- read.table(file.path(INF,"work","INF1.merge-rsid.cs"),col.names=c("prot","MarkerName","CredibleSet"),sep="\t")
pqtls <- merge(pqtls,credibleset,by.x=c("prot","rsid"),by.y=c("prot","MarkerName")) %>%
         rename(Protein=prot) %>% mutate(prots=Protein,Protein=target.short) %>% select(-target.short)
metal <- read.delim(file.path(INF,"work","INF1.METAL"))
pqtldisease <- subset(read.sheet("short",1:51,1:220),Keep==1) %>%
               left_join(unique(pqtls[c("Protein","prots")]),by="prots") %>%
               mutate(prots=if_else(grepl(";",prots),prots,str_replace(prots,prots,Protein))) %>%
               mutate(prots=if_else(grepl("CXCL9;IL.12B",prots),str_replace_all(prots,c("IL.12B"="IL-12B")),prots)) %>%
               mutate(prots=if_else(grepl("MMP.10;CST5",prots),str_replace_all(prots,c("MMP.10"="MMP-10")),prots)) %>%
               rename(Proteins=prots) %>% left_join(distinct(metal[c("MarkerName","rsid")])) %>%
               select(rsid,Proteins,Allele1,Allele2,Effects,SEs,cistrans,trait,efo,study,pmid,dataset,infection)
coloc <- merge(read.delim(file.path(INF,"coloc","GTEx-all.tsv")),gap_inf1,by="prot") %>%
         mutate(prot=target.short,flag=if_else(H3+H4>=0.9 & H4/H3>=3,"x","")) %>%
         rename(Protein=prot) %>% select(-target.short) %>% arrange(desc(flag))
cs95 <- read.delim(file.path(INF,"coloc.1M","cis-eQTL_table.tsv"))
cs95 <- data.frame(rsidProt=str_replace(rownames(cs95),"[.]","-"),cs95)
HOME <- Sys.getenv("HOME")
load(file.path(HOME,"software-notes","docs","files","pi_database.rda"))
drug <- subset(pi_drug,target%in%with(gap.datasets::inf1,gene)) %>% left_join(pi_trait)
efo <- read.delim(file.path(INF,"rsid","efo.txt"))
hgi <- read.delim(file.path(INF,"mr","gsmr","hgi","5e-8","5e-8.tsv"))
pqtls <- select(pqtls,-prots)

outsheets <- c("summary","studies","inf1",
               "pqtls","cojo","knownpqtls","coloc","cs95","pqtldisease",
               "vep","garfield",
               "gsmr_efo","hgi","drug",
               "reactome","great","efo","gdb",
               "interval","os","cvd1","aristotl","pqtlstudies",
               "great3","mr_immun","smr","mr","mr_misc","gsmr",
               "protein_correlation", "protein_dgi", "pqtl_impact")
titles <- c("summary","study information","panel information",
            "pQTLs","conditional analysis","known pQTLs","GTEx coloc","GTEx coloc 95%CS","Disease GWAS overlap",
            "VEP annotation","GARFIELD outputs",
            "GSMR results","HGI r6","PI drug",
            "Reactome","GREAT","EFO","geneDrugbank",
            "INTERVAL study","Other studies","SCALLOP-CVD1","ARISTOTLE study","previous pQTL studies",
            "IL12B-KITLG-TNFSF10","pQTL-immune-MR","SMR","MR results","pQTL-misc-MR","GSMR-FEV1CVD",
            "Protein correlation","DGI membership", "pQTL impact")
description=paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))
uppered <- c("PQTLs")
description[description%in%uppered] <- titles[description%in%uppered]
n0 <- 3
n1 <- 11
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
# if (i<=n0+n1 | i>n0+n1+3)
# {
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
# } else {
#   sheet <- paste0("Annotrans-",i-n0-n1)
#   addWorksheet(wb,sheet,gridLines=FALSE)
#   writeData(wb,sheet,paste0("Trans-pQTL annotation-",i-n0-n1), xy=c(1,1),
#             headerStyle=createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
#   writeDataTable(wb,sheet, get(paste0("at",i-n0-n1)), xy=c(1,2), firstColumn=TRUE, bandedRows=FALSE)
# }
}
sheets_wb <- sheets(wb)
data.frame(sheets_wb)

bStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
hStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
#conditionalFormatting(wb, sheets_wb[grepl("IL12B-KITLG-TNFSF10",sheets_wb)], cols = 26, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GREAT",sheets_wb)], cols = 25, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GTEx coloc$",sheets_wb)], cols = 12, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GARFIELD",sheets_wb)], cols = 3, rows = 3:nrow(garfield), rule = "<=1e-5", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("immune-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_immun), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("MR results",sheets_wb)], cols = 9, rows = 3:nrow(mr), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("misc-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_misc), rule = "==\"x\"", style = hStyle)

saveWorkbook(wb, file=xlsx, overwrite=TRUE)

# mr_immun <- read.sheet("pqtlMR-immune", 1:7, 2:67)
#  mr_misc <- read.sheet("pqtlMR-misc", 1:7, 2:39)
#      ivw <- read.sheet("IVW", 1:8, 2:19)
#     mrc2 <- read.sheet("MRC2", 1:7, 2:8)
#     mvmr <- read.sheet("MVMR", 1:9, 2:6)

novelpqtls <- subset(within(pqtls,{
                                    chrpos=paste0(Chromosome,":",Position)
                                    a1a2=paste0(toupper(Allele1),"/",toupper(Allele2))
                                    bse=paste0(round(Effect,3)," (",round(StdErr,3),")")
                                    log10p=-log.P.
                                  }),
                     !paste0(Protein,"-",rsid)%in%with(knownpqtls,paste0(Protein,"-",Sentinels)),
                     select=c(Protein,uniprot,chrpos,rsid,a1a2,bse,log10p,cis.trans)) %>%
              arrange(Protein,chrpos)
write.xlsx(cbind(no=1:nrow(novelpqtls),novelpqtls), file=file.path(INF,"NG","novelpqtls.xlsx"), overwrite=TRUE,
           colNames=TRUE,
           borders="surrounding", headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")

suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
# Input
options(width=200)
HPC_WORK <- Sys.getenv("HPC_WORK")
INF <- Sys.getenv("INF")
rt <- file.path(INF,"mr","gsmr")

xlsx <- "https://jhz22.user.srcf.net/INF/latest/efo-update.xlsx"
efo_update <- read.xlsx(xlsx,sheet=1) %>%
              filter(!is.na(Disease)) %>%
              mutate(opengwasid=unlist(lapply(strsplit(Source,"/"),"[",5))) %>%
              arrange(opengwasid)
efo_info <- ieugwasr::gwasinfo(pull(efo_update,opengwasid)) %>%
            mutate(trait=gsub("\\b(^[a-z])","\\U\\1",substr(trait,1,30),perl=TRUE)) %>%
            select(id,trait,unit,pmid,author,year,ncase,ncontrol,nsnp,population)
knitr::kable(efo_info)

csv <- read.delim(file.path(INF,"OpenGWAS","finngen_endpoints.tsv"))
finngen <- with(efo_update,grepl("finn-b",opengwasid))
finngen_r7_N <- subset(csv,phenocode %in%sel.trait)[sel.var] %>%
                arrange(phenocode) %>%
                select(number.of.cases, number.of.controls) %>%
                mutate(all=number.of.cases+number.of.controls)
efo_update[finngen,c("N.cases","N.controls","N.total")] <- finngen_r7_N
write.table(efo_update,file=file.path(INF,"OpenGWAS","efo-update.txt"),quote=FALSE,row.names=FALSE,sep="\t")
knitr::kable(efo_update)
knitr::kable(cbind(efo_update[c(1:3,6)],efo_info[c(2:4,7:9)]))
sel.trait <- c("D3_SARCOIDOSIS","L12_PSORIASIS","M13_ANKYLOSPON")
sel.var <- c("phenotype","phenocode","number.of.cases","number.of.controls")
knitr::kable(subset(csv,phenocode %in%sel.trait)[sel.var] %>% arrange(phenotype))

efo_2022_10_05 <- function()
{
  xlsx <- "https://jhz22.user.srcf.net/INF/latest/GSMR_datasets.xlsx"
  efo <- read.xlsx(xlsx,sheet=1,startRow=2,colNames=TRUE,skipEmptyRows=TRUE)
  efo_left <- efo[1:5] %>%
              rename(Source1=Source) %>%
              mutate(N.cases=as.numeric(N.cases),N.controls=as.numeric(N.controls),Total.N=as.numeric(Total.N))
  efo_right <- efo[6:10] %>% rename(Source2=Source) %>% mutate(N.cases=as.numeric(N.cases),N.controls=as.numeric(N.controls),Total.N=as.numeric(Total.N))
  efo_old <- filter(cbind(efo_left,efo_right[5]),!grepl("http",Source2)) %>%
             select(-Source2) %>% rename(Source=Source1)
  efo_new <- filter(efo_right,grepl("http",Source2)) %>%
             rename(Source=Source2)
  efo_update <- bind_rows(efo_old,efo_new) %>%
                filter(!is.na(Disease))
  missed <- with(efo_update, grepl("Crohn",Disease))
  efo_update[missed,c("N.cases","N.controls","Source")] <- c(12194,280722,"https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004132/")
  ra <- with(efo_update, grepl("Rheumatoid",Disease))
  efo_update[ra,c("N.cases","N.controls","Source")] <- c(19234,61565,"https://gwas.mrcieu.ac.uk/datasets/ieu-a-833")
  efo_update <- mutate(efo_update,opengwasid=unlist(lapply(strsplit(Source,"/"),"[",5)))
  efo_update <- arrange(efo_update,opengwasid)
  sel.trait <- c("JUVEN_ARTHR","D3_SARCOIDOSIS","L12_PSORIASIS","M13_SLE")
  finngen <- with(efo_update,grepl("finn-b",opengwasid))
  efo_update[finngen,c("N.cases","N.controls","Total.N")]
  finngen_r7_N <- subset(csv,phenocode %in%sel.trait)[sel.var] %>%
                  arrange(phenocode) %>%
                  select(number.of.cases, number.of.controls) %>%
                  mutate(all=number.of.cases+number.of.controls)
  efo_update[finngen,c("N.cases","N.controls","Total.N")] <- finngen_r7_N
  write.table(efo_update,file="efo-update.txt",quote=FALSE,row.names=FALSE,sep="\t")
  knitr::kable(efo_update)
  # Output
  xlsx <- file.path("efo-update.xlsx")
  wb <- createWorkbook(xlsx)
  hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
  for (sheet in c("efo_update"))
  {
     addWorksheet(wb,sheet,zoom=150)
     writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
               fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
     body <- get(sheet)
     writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
     freezePane(wb, sheet, firstCol=TRUE, firstActiveRow=3)
     width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
   # width_vec_header <- nchar(colnames(body))+2
     setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
     writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
  }
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}

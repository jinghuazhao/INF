suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
# Input
options(width=200)
HPC_WORK <- Sys.getenv("HPC_WORK")
INF <- Sys.getenv("INF")
rt <- file.path(INF,"mr","gsmr")
prot <- "TNFB"
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
efo_update[missed,c("N.cases","N.controls","Source")] <- c("12194","28072","https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004132/")
efo_update <- mutate(efo_update,opengwasid=unlist(lapply(strsplit(Source,"/"),"[",5)))
write.table(efo_update,file="efo_update.txt",quote=FALSE,row.names=FALSE,sep="\t")
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

library(TwoSampleMR)
library(ggplot2)
library(stringr)

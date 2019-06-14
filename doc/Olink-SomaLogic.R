# 14-6-2019 JHZ

library(reshape)
Olink <- "INF/doc/olink.inf.panel.annot.tsv"
SomaLogic <- "SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv"

o <- read.delim(paste0("/home/jhz22/",Olink))
s <- read.delim(paste0("/home/jhz22/",SomaLogic), as.is=TRUE)
s <- rename(s, c(UniProt="uniprot"))

os <- merge(o,s,by="uniprot")

library(dplyr)
nj <- nest_join(o,s,by="uniprot")
unique(unlist(lapply(nj[["y"]],"[",7)))

xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:19), rows=c(3:2746))
t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:20), rows=c(3:167))


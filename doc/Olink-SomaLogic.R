# 13-6-2019 JHZ

library(reshape)
Olink <- "INF/doc/olink.inf.panel.annot.tsv"
SomaLogic <- "SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv"

o <- read.delim(paste0("/home/jhz22/",Olink))
s <- read.delim(paste0("/home/jhz22/",SomaLogic), as.is=TRUE)
s <- rename(s, c(UniProt="uniprot"))

os <- merge(o,s,by="uniprot")

library(dplyr)
nj <- nest_join(o,s,by="uniprot")

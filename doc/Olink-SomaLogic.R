# 9-7-2020 JHZ

Olink <- "INF/doc/olink.inf.panel.annot.tsv"
o <- read.delim(paste0("/home/jhz22/",Olink))
SomaLogic <- "SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv"
s <- read.delim(paste0("/home/jhz22/",SomaLogic), as.is=TRUE)

library(reshape)
s <- rename(s, c(UniProt="uniprot"))
setdiff(intersect(o["uniprot"],s["uniprot"]),"P23560")

os <- merge(o,s,by="uniprot")
u <- setdiff(unique(os$uniprot),"P23560")
write(u,file="u")
length(u)
unique(subset(os[c("uniprot","target.short")],uniprot%in%u))

library(dplyr)
nj <- nest_join(o,s,by="uniprot")
nj[order(nj["uniprot"]),"uniprot"]

library(VennDiagram)
plist <- list(setdiff(o[["uniprot"]],"P23560"),setdiff(s[["uniprot"]],c(NA,"P23560")))
cnames=c("Olink", "SomaLogic")
venn.diagram(x = plist, category.names=cnames, filename='Olink-SomaLogic.png', imagetype="png", output=TRUE)

xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:19), rows=c(3:2746))
t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:20), rows=c(3:167))

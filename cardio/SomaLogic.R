# 30-4-2019 JHZ

SomaLogic <- read.delim("/home/jhz22/SomaLogic/doc/SOMALOGIC_Master_Table_160410_1129info.tsv",as.is=TRUE)
INF <- read.delim("/home/jhz22/INF/doc/olink.inf.panel.annot.tsv",as.is=TRUE)
SomaLogic_INF <- merge(SomaLogic,INF,by.x="UniProt",by.y="uniprot")
t <- with(SomaLogic_INF,table(target.short))
length(t)
names(t)

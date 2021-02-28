# Olink/SomaLogic overlap

library(pQTLtools)
library(VennDiagram)
library(dplyr)
HOME <- Sys.getenv("HOME")
INF <- Sys.getenv("INF")

# All the in the SomaLogic/Olink panels
head(SomaLogic160410)
head(inf1)
A <- setdiff(with(SomaLogic160410,UniProt),"P23560")
B <- setdiff(with(inf1,uniprot),"P23560")
# 3523, 77, 14
venn.diagram(list(SomaLogic=A, Olink=B),"uniprot.tiff",height=8,width=8,units="in")

# Those in common and with pQTLs
overlap <- intersect(A,B)
INF1 <- read.delim(file.path(INF,"work","INF1.METAL")) %>% filter(uniprot %in% overlap) %>% select(uniprot)
# pQTLs from Box
# 15, 54, 5
INTERVAL_box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv")) %>% filter(UniProt %in% overlap) %>% select(UniProt)
venn.diagram(list(SomaLogic=unique(INTERVAL[["UniProt"]]), Olink=unique(INF1[["uniprot"]])),"box-INF1.tiff",height=8,width=8,units="in")
# pQTLs from PhenoScanner
# 4, 23, 36
INTERVAL <- read.delim(file.path(INF,"ps","pQTL.Sun-B_pQTL_EUR_2017.tsv")) %>% filter(UniProt %in% overlap) %>% select(UniProt)
venn.diagram(list(SomaLogic=unique(INTERVAL[["UniProt"]]), Olink=unique(INF1[["uniprot"]])),"ps-INF1.tiff",height=8,width=8,units="in")
# ST4
# 5, 27, 32
ST4 <- st4 %>% select(UniProt) %>% filter(UniProt %in% overlap)
venn.diagram(list(SomaLogic = unique(ST4[["UniProt"]]), Olink = unique(INF1[["uniprot"]])),"ST4-INF1.tiff",height=8,width=8,units="in")
cistrans <- st4 %>% select("cis/.trans",UniProt) %>% rename(cistrans="cis/.trans") %>% filter(UniProt %in% overlap)
# 16, 43
cis <- cistrans %>% filter(cistrans == "cis") %>% select(UniProt)
B <- INF1_metal %>% filter(cis.trans=="cis") %>% select(uniprot) %>% unique
venn.diagram(list(SomaLogic = unique(cis[["UniProt"]]), Olink = B[["uniprot"]]),"cis.tiff",height=8,width=8,units="in")
trans <- cistrans %>% filter(cistrans == "trans") %>% select(UniProt)
B <- INF1_metal %>% filter(cis.trans=="trans") %>% select(uniprot) %>% unique
# 8, 11, 40
venn.diagram(list(SomaLogic = unique(trans[["UniProt"]]), Olink = B[["uniprot"]]),"trans.tiff",height=8,width=8,units="in")
# coloc???

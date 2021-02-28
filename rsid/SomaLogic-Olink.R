# Olink/SomaLogic overlap

library(pQTLtools)
library(VennDiagram)
library(dplyr)
HOME <- Sys.getenv("HOME")
INF <- Sys.getenv("INF")

source.vs.INF1 <- function(source,group)
  venn.diagram(list(SomaLogic=unique(source[["UniProt"]]), Olink=INF1[["uniprot"]]),paste0(group,"-INF1.tiff"),
               force.unique=FALSE,height=8,width=8,units="in")
# panels
SomaLogic <- setdiff(with(SomaLogic160410,UniProt),"P23560")
Olink <- setdiff(with(inf1,uniprot),"P23560")
venn.diagram(list(SomaLogic=SomaLogic, Olnk=Olink),"SomaLogic-Olink.tiff",force.unique=FALSE,height=8,width=8,units="in")
overlap <- intersect(SomaLogic,Olink)
# pQTLs
INF1_metal <- read.delim(file.path(INF,"work","INF1.METAL"))
INF1 <- INF1_metal %>% filter(uniprot %in% overlap)
box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv")) %>% filter(UniProt %in% overlap)
source.vs.INF1(box,"box")
ps <- read.delim(file.path(INF,"ps","pQTL.Sun-B_pQTL_EUR_2017.tsv")) %>% filter(UniProt %in% overlap)
source.vs.INF1(ps,"ps")
ST4 <- st4 %>% select(Locus.ID,SOMAmer.ID,Target,UniProt,"Sentinel.variant*",Chr,Pos,"cis/.trans") %>%
       filter(UniProt %in% overlap) %>%
       rename(cis_trans="cis/.trans") %>% mutate(UniProt_group=paste0(UniProt,"_",cis_trans))
source.vs.INF1(ST4,"ST4")

INF1 <- INF1_metal %>% filter(uniprot %in% overlap) %>% mutate(uniprot_group=paste0(uniprot,"_",cis.trans))
venn.diagram(list(Olink=with(INF1,uniprot),SomaLogic=with(ST4,UniProt)),
             file="ab.tiff",force.unique=FALSE,height=8,width=8,units="in")
venn.diagram(list(Olink=INF1[grepl("cis",INF1$cis.trans),"uniprot"],SomaLogic=ST4[grepl("cis",ST4$cis_trans),"UniProt"]),
             file="cis-cis.tiff",force.unique=FALSE,height=8,width=8,units="in")
venn.diagram(list(Olink=INF1[grepl("trans",INF1$cis.trans),"uniprot"],SomaLogic=ST4[grepl("trans",ST4$cis_trans),"UniProt"]),
             file="trans-trans.tiff",force.unique=FALSE,height=8,width=8,units="in")
venn.diagram(list(Olink=INF1[grepl("cis",INF1$cis.trans),"uniprot"],SomaLogic=ST4[grepl("trans",ST4$cis_trans),"UniProt"]),
             file="cis-trans.tiff",force.unique=FALSE,height=8,width=8,units="in")
venn.diagram(list(Olink=INF1[grepl("trans",INF1$cis.trans),"uniprot"],SomaLogic=ST4[grepl("cis",ST4$cis_trans),"UniProt"]),
             file="trans-cis.tiff",force.unique=FALSE,height=8,width=8,units="in")
# Name	SL	common	Olink
# Panel	3523	77	14
# ---
# Box	15	54	102
# PS	4	23	133
# agree
# cis	5	16	33
# trans	17	11	96
# ST4	22	27	129
# disagree
# c-t	18	10	39
# t-c	8	13	94

# coloc???

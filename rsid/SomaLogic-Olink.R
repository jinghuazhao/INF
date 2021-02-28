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
ST4 <- st4 %>% select(Locus.ID,SOMAmer.ID,Target,UniProt,"Sentinel.variant*",Chr,Pos,"cis/.trans") %>% filter(UniProt %in% overlap) %>%
       rename(cis_trans="cis/.trans") %>% mutate(UniProt_group=paste0(UniProt,"_",cis_trans))
source.vs.INF1(ST4,"ST4")

INF1 <- INF1 %>% mutate(uniprot_group=paste0(uniprot,"_",cis.trans))
cistrans <- ST4 %>% filter(UniProt %in% overlap)
# 22, 27, 129
venn.diagram(list(Olink=with(INF1,uniprot_group),SomaLogic=subset(ST4,UniProt %in% overlap)$UniProt_group),file="ab.tiff",force.unique=FALSE)

cis.vs.trans <- function(A,B,group)
{
  A <- A %>% filter(cis_trans==group) %>% select(UniProt_group)
  B <- B %>% filter(cis.trans==group & uniprot %in% overlap) %>% select(uniprot_group)
  venn.diagram(list(SomaLogic = A[[names(A)]], Olink = B[[names(B)]]),paste0(group,".tiff"),force.unique=FALSE,height=8,width=8,units="in")
}

cis.vs.trans(cistrans,INF1,"cis")
cis.vs.trans(cistrans,INF1,"trans")

ab <-merge(cistrans,INF1,by.x="UniProt",by.y="uniprot")
with(ab,table(cis_trans,cis.trans))
# 99, 27, 99
venn.diagram(as.list(ab[c("UniProt_group","uniprot_group")]),filename="ab-all.tiff",force.unique=FALSE)

# Name	SL	common	Olink
# Panel	3523	77	14
# ---
# Box	15	54	102
# PS	4	23	133
# ---
# cis	5	16	33
# trans	17	11	96
# ST4	22	27	129

#     SL   INF1 (cis.trans)
# cis_trans cis trans  total
#     cis    21    35  56
#     trans  17    53  70
#     total  38    88  126

# coloc ???

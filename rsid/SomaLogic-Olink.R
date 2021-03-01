# Olink/SomaLogic overlap

library(pQTLtools)
library(VennDiagram)
library(dplyr)
HOME <- Sys.getenv("HOME")
INF <- Sys.getenv("INF")

# panels
SomaLogic <- setdiff(with(SomaLogic160410,UniProt),"P23560")
Olink <- setdiff(with(inf1,uniprot),"P23560")
venn.diagram(list(SomaLogic=SomaLogic, Olnk=Olink),"SomaLogic-Olink-panel.tiff",force.unique=FALSE,height=8,width=8,units="in")

# pQTLs
source.vs.INF1 <- function(source,group)
  venn.diagram(list(SomaLogic=source[["UniProt"]], Olink=INF1[["uniprot"]]),paste0(group,"-INF1-pQTLs.tiff"),
               force.unique=FALSE,height=8,width=8,units="in")
overlap <- intersect(SomaLogic,Olink)
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

vd <- function(d,f) venn.diagram(as.list(d),category=c("SomaLogic","Olink"),force.unique=FALSE,height=8,width=8,units="in",filename=f)
join <- function()
{
  INF1 <- INF1_metal %>% filter(uniprot %in% overlap) %>% mutate(uniprot_group=paste0(uniprot,"_",cis.trans))
  ST4_INF1 <-merge(ST4,INF1,by.x="UniProt",by.y="uniprot")
  with(ST4_INF1,table(cis_trans,cis.trans))

  AA <- ST4_INF1[c("uniprot_group","UniProt_group")]
  CC <- subset(ST4_INF1,grepl("cis",cis.trans) & grepl("cis",cis_trans), select=c(cis_trans,cis.trans))
  CT <- subset(ST4_INF1,grepl("cis",cis.trans) & grepl("trans",cis_trans), select=c(cis_trans,cis.trans))
  TC <- subset(ST4_INF1,grepl("trans",cis.trans) & grepl("cis",cis_trans), select=c(cis_trans,cis.trans))
  TT <- subset(ST4_INF1,grepl("trans",cis.trans) & grepl("trans",cis_trans), select=c(cis_trans,cis.trans))

  vd(AA,"ST4-INF1.tiff")
  vd(CC,"cis-cis.tiff")
  vd(CT,"cis-trans.tiff")
  vd(TC,"trans-cis.tiff")
  vd(TT,"trans-trans.tiff")
}

UniProts <- unique(ST4$UniProt)
uniprots <- unique(INF1$uniprot)
both <- intersect(ST4$UniProt,INF1$uniprot)
venn.diagram(list(SomaLogic=ST4[["UniProt"]],Olink=INF1[["uniprot"]]),filename="ST4-INF1-prot.tiff")

CC <- list(subset(ST4,UniProt %in% both & cis_trans=="cis")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="cis")[["uniprot"]])
CT <- list(subset(ST4,UniProt %in% both & cis_trans=="cis")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="trans")[["uniprot"]])
TC <- list(subset(ST4,UniProt %in% both & cis_trans=="trans")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="cis")[["uniprot"]])
TT <- list(subset(ST4,UniProt %in% both & cis_trans=="trans")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="trans")[["uniprot"]])
vd(CC,"cis-cis.tiff")
vd(CT,"cis-trans.tiff")
vd(TC,"trans-cis.tiff")
vd(TT,"trans-trans.tiff")

#
# Name	SL	both	Olink
# proteins
# Panel	3523	77	14
# ST4	5	27	32
# pQTLs
# Box	48	54	102
# PS	312	23	133
# ST4	22	27	129
# ST4/INF1
# CC	5	16	7
# TT	11	11	41
# CT	8	13	39
# TC	12	10	13
#
# coloc???

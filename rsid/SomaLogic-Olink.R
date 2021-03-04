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
       mutate(UniProt=if_else(UniProt=="P29460,Q9NPF7","P29460",UniProt)) %>%
       filter(UniProt %in% overlap) %>%
       rename(cis_trans="cis/.trans")
source.vs.INF1(ST4,"ST4")

INF1 <- INF1_metal %>% filter(uniprot %in% overlap)
UniProts <- unique(ST4$UniProt)
uniprots <- unique(INF1$uniprot)
both <- intersect(ST4$UniProt,INF1$uniprot)
UniProt_uniprot <- list(SomaLogic=ST4[["UniProt"]],Olink=INF1[["uniprot"]])
venn.diagram(UniProt_uniprot,filename="ST4-INF1-prot.tiff",fill=c("blue","red"),height=6,width=6,units="in")
calc.overlap <- calculate.overlap(UniProt_uniprot)
cat(unlist(lapply(calc.overlap,length)),sep="\t","\n")

C1 <- subset(ST4,UniProt %in% both & cis_trans=="cis")[["UniProt"]]
C2 <- subset(INF1,uniprot %in% both & cis.trans=="cis")[["uniprot"]]
T1 <- subset(ST4,UniProt %in% both & cis_trans=="trans")[["UniProt"]]
T2 <- subset(INF1,uniprot %in% both & cis.trans=="trans")[["uniprot"]]
C1C2 <- intersect(C1,C2)
T1T2 <- intersect(T1,T2)

CC <- list(subset(ST4,UniProt %in% both & cis_trans=="cis")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="cis")[["uniprot"]])
TT <- list(subset(ST4,UniProt %in% both & cis_trans=="trans")[["UniProt"]],subset(INF1,uniprot %in% both & cis.trans=="trans")[["uniprot"]])
CT <- list(setdiff(C1,C1C2),setdiff(T2,T1T2))
TC <- list(setdiff(T1,T1T2),setdiff(C2,C1C2))

vd <- function(d,f,diagram=FALSE)
{
  dlist <- as.list(d)
  if (diagram) venn.diagram(dlist,category=c("SomaLogic","Olink"),force.unique=FALSE,height=8,width=8,units="in",filename=f)
  calc.overlap <- calculate.overlap(dlist)
  a1_a2_a3 <- unlist(lapply(calc.overlap,length))
  cat(c(a1_a2_a3[1]-a1_a2_a3[3],a1_a2_a3[3],a1_a2_a3[2]-a1_a2_a3[3]),sep="\t","\n")
}

vd(CC,"cis-cis.tiff", diagram=TRUE)
vd(TT,"trans-trans.tiff", diagram=TRUE)
vd(CT,"cis-trans.tiff")
vd(TC,"trans-cis.tiff")
OnlyUniProts <- setdiff(with(ST4,UniProt),both)
Onlyuniprots <- setdiff(with(INF1,uniprot),both)
table(subset(ST4,UniProt %in% OnlyUniProts)["cis_trans"])
table(subset(ST4,UniProt %in% both)["cis_trans"])
table(subset(INF1,uniprot %in% Onlyuniprots)["cis.trans"])
table(subset(INF1,uniprot %in% both)["cis.trans"])

#
# Name	SL	both	Olink
# proteins
# Panel	3523	77	14
# ST4	5	28	31
# pQTLs
# Box	48	54	102
# PS	312	23	133
# ST4	23	28	128
# ST4/INF1
# CC	5	17	7
# TT	11	12	47
# CT	0	0	10
# TC	3	0	7
#
# coloc???

SomaLogic_Olink <- function()
{
  library(VennDiagram)
  INF1_prot <- read.table("work/INF1.merge.prot",col.names=c("prot","uniprot"))
  library(pQTLtools)
  SomaLogic_prot <- unique(replace(st4$UniProt,st4$UniProt=="P29460,Q9NPF7","P29460"))
  plist <- list(INF1_prot$uniprot,SomaLogic_prot)
  ov <- VennDiagram::calculate.overlap(plist)
  ov$a3
  venn.plot <- draw.pairwise.venn(1469, 70, 28,
               category = c("SomaLogic", "Olink"),
               fill = c("blue", "red"),
               lty = "blank",
               cex = 2,
               cat.cex = 2,
               cat.pos = c(200, 50),
               cat.dist = 0.09,
               cat.just = list(c(-1, -1), c(1, 1)),
               ext.pos = 30,
               ext.dist = -0.05,
               ext.length = 0.85,
               ext.line.lwd = 2,
               ext.line.lty = "dashed"
             );
  grid.draw(venn.plot);
  png('SomaLogic-Olink-proteins.png', height=20, width=20, units="cm", res=300)
  grid.draw(venn.plot);
  dev.off();
  grid.newpage()
  venn.plot <- draw.pairwise.venn(45, 83, 10,
               category = c("SomaLogic", "Olink"),
               fill = c("blue", "red"),
               lty = "blank",
               cex = 2,
               cat.cex = 2,
               cat.pos = c(30, 10),
               cat.dist = 0.09,
               cat.just = list(c(-1, -1), c(1, 1)),
               ext.pos = 30,
               ext.dist = -0.05,
               ext.length = 0.85,
               ext.line.lwd = 2,
               ext.line.lty = "dashed"
             );
  grid.draw(venn.plot);
  png('SomaLogic-Olink-sentinels.png', height=20, width=20, units="cm", res=300)
  grid.draw(venn.plot);
  dev.off();
  grid.newpage()
  st6[c(39,53),"UniProt"] <- "P29460"
  st6ov <- subset(st6,UniProt%in%ov$a3)
# significant on INTERVAL (27)
  dim(subset(st6ov,as.numeric(p.1)<1e-11))
  z <- with(st6ov,{
    chr <- st6ov[["Chr"]]
    pos <- st6ov[["Pos"]]
    a1 <- st6ov[["Effect.Allele.(EA)"]]
    a2 <- st6ov[["Other.Allele.(OA)"]]
    cbind(UniProt,snpid=gap::chr_pos_a1_a2(chr,pos,a1,a2))
  })
  write.table(merge(inf1,z,by.x="uniprot",by.y="UniProt")[c("snpid","prot","uniprot")],
              file="SomaLogic.id3",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

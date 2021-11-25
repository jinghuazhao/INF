# Olink/Somascan overlap

options(width=200)
HOME <- Sys.getenv("HOME")
INF <- Sys.getenv("INF")

library(pQTLtools)
library(VennDiagram)
library(dplyr)

header <- function(title) {plot(1:10, type = "n", ann = FALSE, axes = FALSE);text((10-length(title))/2,5,title)}

# 1. panels
Somascan <- with(SomaLogic160410,unique(UniProt))
Olink <- setdiff(with(inf1,uniprot),"P23560")
overlap <- intersect(Somascan,Olink)
INF1_metal <- read.delim(file.path(INF,"work","INF1.METAL"))
INF1 <- INF1_metal %>% filter(uniprot %in% overlap)
pdf("Somascan-Olink.pdf")
panel <- venn.diagram(list(Somascan=Somascan, Olink=Olink),
                      filename=NULL,force.unique=FALSE,height=8,width=8,
                      main.cex=3,sub.cex=3,cat.cex=3,cat.pos = c(-20, 0),cex=3,units="in")
header("1. Somascan-Olink overlap according to unique uniprot IDs")
grid.newpage()
grid.draw(panel)
# 2. ST4 all proteins
Somascan_prot <- unique(replace(st4$UniProt,st4$UniProt=="P29460,Q9NPF7","P29460"))
INF1_prot <- read.table(file.path(INF,"work","INF1.merge.prot"),col.names=c("prot","uniprot"))
plist <- list(Somascan_prot,INF1_prot$uniprot)
calc.overlap <- calculate.overlap(plist)
a1_a2_a3 <- unlist(lapply(calc.overlap,length))
a1_a2_a3
cat(c(a1_a2_a3[1]-a1_a2_a3[3],a1_a2_a3[3],a1_a2_a3[2]-a1_a2_a3[3]),sep="\t");cat("\n")
header("2. ST4 proteins regardless overlap")
grid.newpage()
venn.plot <- draw.pairwise.venn(a1_a2_a3[1],a1_a2_a3[2],a1_a2_a3[3],
             category = c("Somascan", "Olink"),
             fill = c("blue", "red"),
             lty = "blank",
             cex = 3,
             cat.cex = 3,
             cat.pos = c(-20, 0)
           );
grid.draw(venn.plot)
# 3. ST4 overlapping proteins, Locus.ID=6_10 should be merged
ST4 <- st4 %>% select(Locus.ID,SOMAmer.ID,Target,UniProt,"Sentinel.variant*",Chr,Pos,"cis/.trans") %>%
       mutate(UniProt=if_else(UniProt=="P29460,Q9NPF7","P29460",UniProt)) %>%
       filter(UniProt %in% overlap & SOMAmer.ID!="VEGFA.2597.8.3") %>%
       rename(cis_trans="cis/.trans",rsid="Sentinel.variant*")
both <- intersect(ST4[["UniProt"]],INF1[["uniprot"]])
UniProt_uniprot <- list(Somascan=ST4[["UniProt"]],Olink=INF1[["uniprot"]])
ST4_INF1_prot <- venn.diagram(UniProt_uniprot,filename=NULL,force.unique=TRUE,
                              cat.cex=3.5,cat.pos=c(-20,0),cex=3.5,fill=c("blue","red"),height=6,width=6,units="in")
header("3. ST4 Proteins")
grid.newpage()
grid.draw(ST4_INF1_prot)

# pQTLs
source.vs.INF1 <- function(source,group)
{
  dlist <- list(Somascan=source[["UniProt"]], Olink=INF1[["uniprot"]])
  calc.overlap <- calculate.overlap(dlist)
  a1_a2_a3 <- unlist(lapply(calc.overlap,length))
  cat(c(a1_a2_a3[1]-a1_a2_a3[3],a1_a2_a3[3],a1_a2_a3[2]-a1_a2_a3[3]),sep="\t","\n")
  pqtl <- venn.diagram(dlist,filename=NULL,
                       cat.cex=3,cat.pos=c(-20,0),cex=3,force.unique=FALSE,height=8,width=8,units="in")
  grid.newpage()
  grid.draw(pqtl)
}
# 4.
header("4. pQTLs from ST4 overlapped proteins")
source.vs.INF1(ST4,"ST4")
# 5.
header("5. pQTLs from box")
box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv")) %>% filter(UniProt %in% overlap)
source.vs.INF1(box,"box")
# 6.
header("6. pQTLs from PhenoScanner")
ps <- read.delim(file.path(INF,"ps","pQTL.Sun-B_pQTL_EUR_2017.tsv")) %>% filter(UniProt %in% overlap)
source.vs.INF1(ps,"ps")
## ST4
### Total
OnlyUniProts <- setdiff(with(ST4,UniProt),both)
Onlyuniprots <- setdiff(with(INF1,uniprot),both)
table(subset(ST4,UniProt %in% OnlyUniProts)["cis_trans"])
table(subset(ST4,UniProt %in% both)["cis_trans"])
table(subset(INF1,uniprot %in% both)["cis.trans"])
table(subset(INF1,uniprot %in% Onlyuniprots)["cis.trans"])
### Both
ST4_unique <- ST4 %>% group_by(UniProt) %>%
              mutate(c_t=if_else(cis_trans=="cis","C","T")) %>%
              summarize(ST4_cistrans=paste(paste0(c_t,"-",rsid),collapse=",")) %>%
              select(UniProt,ST4_cistrans) %>%
              filter(UniProt %in% both) %>%
              data.frame
INF1_unique <- INF1 %>% group_by(uniprot) %>%
               mutate(c.t=if_else(cis.trans=="cis","C","T")) %>%
               summarize(INF1_cistrans=paste(paste0(c.t,"-",rsid),collapse=",")) %>%
               select(uniprot,INF1_cistrans) %>%
               filter(uniprot %in% both) %>%
               data.frame
ST4_INF1 <- merge(ST4_unique,INF1_unique,by.x="UniProt",by.y="uniprot") %>% data.frame
ST4_INF1_A <- ST4_INF1 %>%
              group_by(UniProt) %>%
              summarize(intersect=paste(intersect(unlist(strsplit(ST4_cistrans,",")),unlist(strsplit(INF1_cistrans,","))),collapse=",")) %>%
              right_join(ST4_INF1) %>%
              data.frame
write.csv(ST4_INF1_A,file="ST4_INF1.csv",row.names=FALSE)

C1 <- subset(ST4,UniProt %in% both & cis_trans=="cis")
T1 <- subset(ST4,UniProt %in% both & cis_trans=="trans")
C2 <- subset(INF1,uniprot %in% both & cis.trans=="cis")
T2 <- subset(INF1,uniprot %in% both & cis.trans=="trans")
C1C2 <- intersect(C1[["UniProt"]],C2[["uniprot"]])
T1T2 <- intersect(T1[["UniProt"]],T2[["uniprot"]])

CC <- list(C1[["UniProt"]],C2[["uniprot"]])
TT <- list(T1[["UniProt"]],T2[["uniprot"]])
CT <- list(setdiff(C1[["UniProt"]],C1C2),setdiff(T2[["uniprot"]],T1T2))
TC <- list(setdiff(T1[["UniProt"]],T1T2),setdiff(C2[["uniprot"]],C1C2))

vd <- function(dlist,diagram=FALSE)
{
  calc.overlap <- calculate.overlap(dlist)
  a1_a2_a3 <- unlist(lapply(calc.overlap,length))
  cat(c(a1_a2_a3[1]-a1_a2_a3[3],a1_a2_a3[3],a1_a2_a3[2]-a1_a2_a3[3]),sep="\t");cat("\n")
  if (diagram)
  {
    vd <- venn.diagram(dlist,category=c("Somascan","Olink"),
                       filename=NULL, force.unique=FALSE,height=8,width=8,units="in",main.cex=2,sub.cex=2,cat.cex=2,cat.pos = c(-20, 0),cex=2)
    grid.newpage()
    grid.draw(vd)
  }
}

# 7.
header("7. cis pQTLs from ST4")
vd(CC, diagram=TRUE)
# 8.
header("8. trans pQTLs from ST4")
vd(TT, diagram=TRUE)
vd(CT)
vd(TC)
# 9.
header("9. ST4 pQTLs")
grid.newpage()
venn.plot <- draw.pairwise.venn(44, 83, 9,
             category = c("Somascan", "Olink"),
             fill = c("blue", "red"),
             lty = "blank",
             cex = 3,
             cat.cex = 3,
             cat.pos = c(0, -20)
           );
grid.draw(venn.plot);
dev.off()

#
# Name	SS	both	OL
# Proteins
# Panel	3523	77	14
# ST4
# All 	1441    28      42
# +pQTL	5	28	31
# pQTLs
# Box	48	54	102
# PS	312	23	133
# ST4	22	28	128
# ST4/INF1
# CC	5	17	7
# TT	11	12	47
# CT	0	0	10
# TC	3	0	7
#

# INTERVAL Somascan data
# ~/rds/post_qc_data/interval/phenotype/somalogic_proteomics

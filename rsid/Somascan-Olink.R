# Olink/Somascan overlap

options(width=200)
HOME <- Sys.getenv("HOME")
INF <- Sys.getenv("INF")

library(pQTLtools)
library(VennDiagram)
library(dplyr)

# panels
# 1.
Somascan <- with(SomaLogic160410,unique(UniProt))
Olink <- setdiff(with(inf1,uniprot),"P23560")
overlap <- intersect(Somascan,Olink)
INF1_metal <- read.delim(file.path(INF,"work","INF1.METAL"))
INF1 <- INF1_metal %>% filter(uniprot %in% overlap)
pdf("Somascan-Olink.pdf")
panel <- venn.diagram(list(Somascan=Somascan, Olink=Olink),
                      filename=NULL,force.unique=FALSE,height=8,width=8,
                      main.cex=3,sub.cex=3,cat.cex=3,cat.pos = c(-20, 0),cex=3,units="in")
grid.newpage()
grid.draw(panel)
grid.text("Somascan-Olink overlap according to unique uniprot IDs")
# 2.
Somascan_prot <- unique(replace(st4$UniProt,st4$UniProt=="P29460,Q9NPF7","P29460"))
INF1_prot <- read.table(file.path(INF,"work","INF1.merge.prot"),col.names=c("prot","uniprot"))
plist <- list(Somascan_prot,INF1_prot$uniprot)
calc.overlap <- calculate.overlap(plist)
a1_a2_a3 <- unlist(lapply(calc.overlap,length))
a1_a2_a3
cat(c(a1_a2_a3[1]-a1_a2_a3[3],a1_a2_a3[3],a1_a2_a3[2]-a1_a2_a3[3]),sep="\t");cat("\n")
grid.newpage()
venn.plot <- draw.pairwise.venn(a1_a2_a3[1],a1_a2_a3[2],a1_a2_a3[3],
             category = c("Somascan", "Olink"),
             fill = c("blue", "red"),
             lty = "blank",
             cex = 3,
             cat.cex = 3,
             cat.pos = c(-20, 0)
           );
grid.draw(venn.plot);
grid.text("ST4 proteins regardless overlap")
# 3.
ST4 <- st4 %>% select(Locus.ID,SOMAmer.ID,Target,UniProt,"Sentinel.variant*",Chr,Pos,"cis/.trans") %>%
       mutate(UniProt=if_else(UniProt=="P29460,Q9NPF7","P29460",UniProt)) %>%
       filter(UniProt %in% overlap) %>%
       rename(cis_trans="cis/.trans",rsid="Sentinel.variant*")
both <- intersect(ST4[["UniProt"]],INF1[["uniprot"]])
UniProt_uniprot <- list(Somascan=ST4[["UniProt"]],Olink=INF1[["uniprot"]])
ST4_INF1_prot <- venn.diagram(UniProt_uniprot,filename=NULL,force.unique=TRUE,
                              cat.cex=3.5,cat.pos=c(-20,0),cex=3.5,fill=c("blue","red"),height=6,width=6,units="in")
grid.newpage()
grid.draw(ST4_INF1_prot)
grid.text("ST4 Proteins")

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
box <- read.delim(file.path(HOME,"SomaLogic","doc","INTERVAL-box.tsv")) %>% filter(UniProt %in% overlap)
# 4.
source.vs.INF1(box,"box")
grid.text("pQTLs from box")
ps <- read.delim(file.path(INF,"ps","pQTL.Sun-B_pQTL_EUR_2017.tsv")) %>% filter(UniProt %in% overlap)
# 5.
source.vs.INF1(ps,"ps")
grid.text("pQTLs from PhenoScanner")
# 6.
source.vs.INF1(ST4,"ST4")
grid.text("pQTLs from ST4")
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
vd(CC, diagram=TRUE);grid.text("cis pQTLs from ST4")
# 8.
vd(TT, diagram=TRUE);grid.text("trans pQTLs from ST4")
vd(CT)
vd(TC)
# 9.
grid.newpage()
venn.plot <- draw.pairwise.venn(45, 83, 10,
             category = c("Somascan", "Olink"),
             fill = c("blue", "red"),
             lty = "blank",
             cex = 3,
             cat.cex = 3,
             cat.pos = c(0, -20)
           );
grid.draw(venn.plot);
grid.text("ST4 pQTLs")

## ST6
st6[c(39,53),"UniProt"] <- "P29460"
st6ov <- subset(st6,UniProt%in%calc.overlap$a3)
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
# ST4	23	28	128
# ST4/INF1
# CC	5	17	7
# TT	11	12	47
# CT	0	0	10
# TC	3	0	7
#

# INTERVAL Somascan data
# ~/rds/post_qc_data/interval/phenotype/somalogic_proteomics

fenland <- function()
{
# Fenland data
  cross_plat <- "https://www.biorxiv.org/content/biorxiv/early/2021/03/19/2021.03.18.435919/DC2/embed/media-2.xlsx?download=true"
  st1 <- openxlsx::read.xlsx(cross_plat, sheet=2, startRow=2, colNames=TRUE)
  somascan_olink <- subset(st1,Olink.panel=="Olink INFLAMMATION(v.3012)")
  st2 <- openxlsx::read.xlsx(cross_plat, sheet=3, startRow=3, colNames=TRUE)
}

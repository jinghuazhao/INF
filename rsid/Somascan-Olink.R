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
# 3. ST4 overlapping proteins
ST4 <- st4 %>% select(Locus.ID,SOMAmer.ID,Target,UniProt,"Sentinel.variant*",Chr,Pos,"cis/.trans") %>%
       mutate(UniProt=if_else(UniProt=="P29460,Q9NPF7","P29460",UniProt)) %>%
       filter(UniProt %in% overlap) %>%
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
venn.plot <- draw.pairwise.venn(45, 83, 10,
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
# ST4	23	28	128
# ST4/INF1
# CC	5	17	7
# TT	11	12	47
# CT	0	0	10
# TC	3	0	7
#

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

# INTERVAL Somascan data
# ~/rds/post_qc_data/interval/phenotype/somalogic_proteomics

# --- protein overlap
# Somalogic proteins with sentinels (1469) - NOTE P29460,Q9NPF7 in SomaLogic
sed '1d' SomaLogic.sentinels | awk 'a[$2]++==0'| wc -l
cut -f2 work/inf1.tmp | grep -v P23560 > work/INF1.uniprot

# number of proteins with sentinels in both Olink and SomaLogic (28)
cut -f2 work/INF1.merge.prot | grep -f - SomaLogic.sentinels | cut -d' ' -f2 | sort | uniq | wc -l
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - work/INF1.merge.prot | wc -l

# --- signal overlap
# all SomaLogic signals in Olink (51)
join -j2 <(sort -k2,2 work/inf1.tmp | grep -v P23560) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) | wc -l
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f work/INF1.uniprot - | wc -l

# all SomaLogic signals from overlapping proteins (45)
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f work/INF1.merge.uniprot - | wc -l
join -j2 <(sort -k2,2 work/INF1.merge.prot) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) > SomaLogic.INF1.all
(
cut -d' ' -f1-3,5,7,8 SomaLogic.INF1.all | \
parallel -C' ' '
  zgrep -w {3} METAL/{2}-1.tbl.gz
  gunzip -c METAL/{2}-1.tbl.gz | \
  awk -vchr={4} -vstart={5} -vend={6} "NR==1||(\$1==chr&&\$2>=start&&\$2<=end&&\$12<-9.30103)" | \
  cut -f1-5,10-12 | \
  gzip -f > INF1.SomaLogic.{1}-{2}-{3}.gz
  export lines=$(gunzip -c INF1.SomaLogic.{1}-{2}-{3}.gz | wc -l | cut -d" " -f1)
  if [ ${lines} -eq 1 ]; then rm INF1.SomaLogic.{1}-{2}-{3}.gz; fi
'
) > INF1.SomaLogic.all

# which are also genomewide significant (32)
awk '$12<-9.30103' INF1.SomaLogic.all | wc -l

# identical signals (10)
cat SomaLogic.INF1.all | \
parallel -C' ' 'awk -v prot={2} -v MarkerName={3} "\$5==prot && \$6==MarkerName" work/INF1.merge'
join <(awk '{print $2"-"$3,$0}' SomaLogic.INF1.all | sort -k1,1) <(awk '{print $5"-"$6,$0}' work/INF1.merge | sort -k1,1)

# Olink overlapping proteins
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - work/INF1.merge.prot | \
cut -f1 | grep -f - work/INF1.merge | \
cut -f5 | sort | uniq | wc -l

# --- SomaLogic --> Olink lookup --- NOTE these were not necessary Olink sentinels
# Similar to ST6 of the SomaLogic paper and pQTLtools/inst/scripts/STs.R
# from the replicates how many were from INF1 (41)
(
cat SomaLogic.id3 | \
parallel -C' ' '
  zgrep -w {1} METAL/{2}-1.tbl.gz
'
) > SomaLogic.INF1
# how many among INF1 variants were genomewide significant (28)
wc -l SomaLogic.INF1
awk '$12<-9.30103' SomaLogic.INF1 | wc -l

(
cat SomaLogic.id3 | \
parallel -C' ' '
  zgrep -w {1} METAL/{2}-1.tbl.gz | \
  awk -v prot={2} -v uniprot={3} "{print prot,uniprot,\$0}"
'
) | \
sort -k5,5 | \
join -a1 -15 -e "NA" - work/INTERVAL.rsid > SomaLogic.INF1-rsid
awk '$14<-9.30103 {print $2, $21}' SomaLogic.INF1-rsid

rm SomaLogic.id3 SomaLogic.INF1 INF1.SomaLogic*gz SomaLogic.INF1.all  SomaLogic.INF1-rsid  SomaLogic.sentinels INF1.SomaLogic.all

epic_fenland <- function()
{
# EPIC-Fenland data
  cross_plat <- "https://www.biorxiv.org/content/biorxiv/early/2021/03/19/2021.03.18.435919/DC2/embed/media-2.xlsx?download=true"
  st1 <- openxlsx::read.xlsx(cross_plat, sheet=2, startRow=2, colNames=TRUE)
  somascan_olink <- subset(st1,Olink.panel=="Olink INFLAMMATION(v.3012)")
  st2 <- openxlsx::read.xlsx(cross_plat, sheet=3, startRow=3, colNames=TRUE)

  export dir=~/rds/results/public/proteomics/EPIC-Fenland
  export M=1000000

  if [ ! -d ${INF}/epic-fenland ]; then mkdir ${INF}/epic-fenland; fi

  (
    awk '{print "prot","uniprot",$0}' ${dir}/header.txt
    cut -f1-5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | \
    sed '1d' | \
    parallel -C' ' --env dir --env M '
      export MarkerName={1}
      export rsid={2}
      export prot={3}
      export chr={4}
      export pos={5}
      export uniprot={6}
      export start=$((${pos}-${M}))
      export end=$((${pos}+${M}))
      export out=${INF}/epic-fenland/${prot}-${MarkerName}.txt
      (
        cat ${dir}/header.txt
        tabix ${dir}/all.grch37.tabix.gz ${chr}:${start}-${end}
      ) > ${out}
      grep -w ${pos} ${out} | \
      awk -v prot=${prot} -v uniprot=${uniprot} "{print prot,uniprot,\$0}"
    '
  ) > ${INF}/epic-fenland/sentinels.txt

  R --no-save <<\ \ END
    library(pQTLtools)
    library(stringr)
    l <- str_split(SomaLogic160410$SOMAMER_ID,"[.]")
    prot <- unlist(lapply(l,"[",1))
    id3 <- unlist(lapply(l,"[",2))
    id4 <- unlist(lapply(l,"[",3))
    id <- paste(id3,id4,sep="_")
    panel <- cbind(SomaLogic160410,prot,id)
    INF <- Sys.getenv("INF")
    sentinels <- read.table(file.path(INF,"epic-fenland","sentinels.txt"),header=TRUE)
    test <- merge(merge(panel,sentinels,by.x="id",by.y="Somamer"),gap::inf1[c("uniprot","prot","target.short")], by="uniprot")
    subset(test[c("MarkerName","chr.x","chr.y","target.short","prot.y","id","Pvalue")],chr.x==chr.y & Pvalue<5e-8)
  END
}

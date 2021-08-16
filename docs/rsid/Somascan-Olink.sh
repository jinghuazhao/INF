#!/usr/bin/bash

R --no-save -q <<END
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

  library(pQTLtools)
# Locus.ID=6_10 should be merged
  sentinels <- subset(st4,SOMAmer.ID!="VEGFA.2597.8.3")[,5:12]
# Somalogic proteins with sentinels (1469)
  length(unique(sort(st4$UniProt)))
  snpid <- gap::chr_pos_a1_a2(sentinels[,3],sentinels[,4],sentinels[,7],sentinels[,8])
  write.table(cbind(snpid,sentinels),file="SomaLogic.sentinels",quote=FALSE,row.names=FALSE)
END

# --- protein overlap
sed '1d' SomaLogic.sentinels | awk 'a[$2]++==0'| wc -l
cut -f2 ${INF}/work/inf1.tmp | grep -v P23560 > ${INF}/work/INF1.uniprot

# number of proteins with sentinels in both Olink and SomaLogic (28) - NOTE P29460,Q9NPF7 in SomaLogic
cut -f2 work/INF1.merge.prot | grep -f - SomaLogic.sentinels | cut -d' ' -f2 | sort | uniq | wc -l
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - ${INF}/work/INF1.merge.prot | wc -l

# --- signal overlap
# all SomaLogic signals in Olink (50)
join -j2 <(sort -k2,2 ${INF}/work/inf1.tmp | grep -v P23560) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) | wc -l
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f ${INF}/work/INF1.uniprot - | wc -l

# all SomaLogic signals from overlapping proteins (45)
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f ${INF}/work/INF1.merge.uniprot - | wc -l
join -j2 <(sort -k2,2 ${INF}/work/INF1.merge.prot) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) > SomaLogic.INF1.all
(
cut -d' ' -f1-3,5,7,8 SomaLogic.INF1.all | \
parallel -C' ' '
  zgrep -w {3} ${INF}/METAL/{2}-1.tbl.gz
  gunzip -c ${INF}/METAL/{2}-1.tbl.gz | \
  awk -vchr={4} -vstart={5} -vend={6} "NR==1||(\$1==chr&&\$2>=start&&\$2<=end&&\$12<-9.30103)" | \
  cut -f1-5,10-12 | \
  gzip -f > INF1.SomaLogic.{1}-{2}-{3}.gz
  export lines=$(gunzip -c INF1.SomaLogic.{1}-{2}-{3}.gz | wc -l | cut -d" " -f1)
  if [ ${lines} -eq 1 ]; then rm INF1.SomaLogic.{1}-{2}-{3}.gz; fi
'
) > INF1.SomaLogic.all

# which are also genomewide significant (31)
awk '$12<-9.30103' INF1.SomaLogic.all | wc -l

# identical signals (9)
cat SomaLogic.INF1.all | \
parallel -C' ' 'awk -v prot={2} -v MarkerName={3} "\$5==prot && \$6==MarkerName" ${INF}/work/INF1.merge'
join <(awk '{print $2"-"$3,$0}' SomaLogic.INF1.all | sort -k1,1) <(awk '{print $5"-"$6,$0}' ${INF}/work/INF1.merge | sort -k1,1) | \
cut -d' ' -f1 | join - <(awk '{print $3"-"$1,$2,$21}' ${INF}/work/INF1.METAL | sort -k1,1)

# Olink overlapping proteins
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - ${INF}/work/INF1.merge.prot | \
cut -f1 | grep -f - ${INF}/work/INF1.merge | \
cut -f5 | sort | uniq | wc -l

# --- SomaLogic --> Olink lookup --- NOTE these were not necessary Olink sentinels
# Similar to ST6 of the SomaLogic paper and pQTLtools/inst/scripts/STs.R
# from the replicates how many were from INF1 (41)
(
  cat SomaLogic.id3 | \
  parallel -C' ' '
    zgrep -w {1} ${INF}/METAL/{2}-1.tbl.gz
  '
) > SomaLogic.INF1
# how many among INF1 variants were genomewide significant (28)
wc -l SomaLogic.INF1
awk '$12<-9.30103' SomaLogic.INF1 | wc -l

(
  cat SomaLogic.id3 | \
  parallel -C' ' '
    zgrep -w {1} ${INF}/METAL/{2}-1.tbl.gz | \
    awk -v prot={2} -v uniprot={3} "{print prot,uniprot,\$0}"
  '
) | \
sort -k5,5 | \
join -a1 -15 -e "NA" - ${INF}/work/INTERVAL.rsid > SomaLogic.INF1-rsid
awk '$14<-9.30103 {print $2, $21}' SomaLogic.INF1-rsid

rm SomaLogic.id3 SomaLogic.INF1 INF1.SomaLogic*gz SomaLogic.INF1.all  SomaLogic.INF1-rsid  SomaLogic.sentinels INF1.SomaLogic.all

function epic_fenland()
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
    test <- merge(merge(panel,sentinels,by.x="id",by.y="Somamer"),gap.datasets::inf1[c("uniprot","prot","target.short")], by="uniprot")
    subset(test[c("MarkerName","chr.x","chr.y","target.short","prot.y","id","Pvalue")],chr.x==chr.y & Pvalue<5e-8)
  END
}

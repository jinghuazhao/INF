#!/usr/bin/bash

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

R --no-save <<END
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

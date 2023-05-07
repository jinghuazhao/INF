#!/usr/bin/bash

# 7-5-2023 JHZ

function nonR()
{
  if [ ! -d ${INF}/aristotle ]; then mkdir ${INF}/aristotle; fi
  (
    gunzip -c ${INF}/sumstats/ARISTOTLE/*txt.gz | \
    head -1 | \
    awk -vOFS="\t" '{print "Protein",$0}'
    cut -f5,6 ${INF}/work/INF1.merge | \
    sed '1d' | \
    tr '\t' ' ' | \
    parallel -j8 -C' ' '
      echo {1}-{2}
      zgrep -w {2} ${INF}/sumstats/ARISTOTLE/ARISTOTLE.{1}.txt.gz;
      awk -vOFS="\t" -vprot={1} "{print prot,\$0}"
    '
  ) > ${INF}/aristotle/INF1.merge.replication.txt

  # snpid --> rsid
  cd ${INF}/aristotle
  for f in INF1.merge.replication.txt
  do
    awk -vOFS="\t" '{if(NF==1) printf $1 OFS; else print}' ${f} > ${f}-rsid
    (
    cat ${INF}/work/INF1.merge.rsid | \
    parallel --dry-run -C' ' "
      export s={1};
      export r={2};
      sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid
    "
    ) | bash
  done
  for p in 5e-10 5e-8 1e-5 3.9e-4 5e-2
  do
    echo ${p}
    awk 'NF>1' ${INF}/aristotle/INF1.merge.replication.txt-rsid | cut -f12 | awk -v p=${p} '$1<p{print $1}' | wc -l
  done
  cd -
  awk -vOFS="\t" '{if(NR>1) {split($1,a,"-");$1=a[1]};print}' ${INF}/aristotle/INF1.merge.replication.txt-rsid | xsel -i
}

Rscript -e '
  options(width=200)
  library(dplyr)
  INF <- Sys.getenv("INF")
  INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))
  aristotle <- read.delim(file.path(INF,"aristotle","INF1.merge.replication.txt-rsid")) %>%
               rename(n=N,Info=RSQ_IMP,Imputed=IMP)
  INF1_aristotle <- mutate(INF1_METAL,Protein=paste0(prot,"-",rsid)) %>%
                    left_join(aristotle) %>%
                    select(Protein,prot,Chromosome,Position,rsid,Allele1,Allele2,EFFECT_ALLELE,REFERENCE_ALLELE,Freq1,Effect,StdErr,log.P.,
                           CODE_ALL_FQ,BETA,SE,PVAL,cis.trans,Info,Imputed) %>%
                    left_join(pQTLdata::inf1[c("prot","target.short")]) %>%
                    mutate(sw=if_else(Allele1==REFERENCE_ALLELE,-1,1),BETA=sw*BETA,sw2=(sign(Effect)==sign(BETA))+0)
  subset(INF1_aristotle[c("Effect","BETA","log.P.","PVAL","cis.trans")],cis.trans=="cis") %>% arrange(Effect)
  table(INF1_aristotle$sw2)
  tbl <- select(INF1_aristotle,target.short,Chromosome,Position,rsid,EFFECT_ALLELE,REFERENCE_ALLELE,CODE_ALL_FQ,
                BETA,SE,PVAL,Info,Imputed,cis.trans) %>%
         rename(Protein=target.short,rsID=rsid,EAF=CODE_ALL_FQ,b=BETA,se=SE,P=PVAL)
  write.table(tbl,file=file.path(INF,"aristotle","aristotle.csv"),row.names=FALSE,quote=FALSE,sep=",")
  png(file=file.path(INF,"aristotle","SF-INF-ARISTOTLE.png"),res=300,width=15,height=15,units="in")
  attach(INF1_aristotle)
  par(mar=c(5,5,1,1))
  plot(Effect,BETA,pch=19,cex=2,col=ifelse(cis.trans=="trans","blue","red"),xaxt="n",ann=FALSE,cex.axis=2)
  axis(1,cex.axis=2,lwd.tick=0.5)
  legend(x=1.5,y=-0.5,c("cis","trans"),box.lwd=0,cex=2,col=c("red","blue"),pch=19)
  mtext("ARISTOTLE",side=2,line=3,cex=2)
  mtext("Meta-analysis",side=1,line=3,cex=2)
  abline(h=0,v=0)
  detach(INF1_aristotle)
  dev.off()
  all <- INF1_aristotle %>%
         filter(!is.na(BETA)) %>%
         select(Effect,BETA)
  cor(all,use="everything")
  cis <- INF1_aristotle %>%
         filter(!is.na(BETA) & cis.trans=="cis") %>%
         select(Effect,BETA)
  cor(cis,use="everything")
  trans <- INF1_aristotle %>%
           filter(!is.na(BETA) & cis.trans=="trans") %>%
           select(Effect,BETA)
  cor(trans,use="everything")
'

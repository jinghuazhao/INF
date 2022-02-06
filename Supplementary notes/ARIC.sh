#!/usr/bin/bash

if [ ! -d ${INF}/ARIC ]; then mkdir ${INF}/ARIC; fi

export sumstats=~/rds/results/public/proteomics/ARIC

(
  cat <(echo -e "seqid\tsnpid\tsnpid38\trsid\tprot") | \
  paste - ${sumstats}/glm.linear.hdr
  cat <(echo pos38) <(cut -f2 ${INF}/work/INF1.b38) | \
  paste ${INF}/work/INF1.METAL - | \
  awk 'NR>1{snpid=$1;gsub(/chr[0-9]*:[0-9]*/,"",$1);print snpid,$0}' | \
  awk '{$2="chr"$5":"$23$2;print snpid,$1,$2,$3,$4,$21}' | \
  sort -k5,5 | \
  join -12 -25 <(sed '1d' ${sumstats}/seqid.txt | cut -f1,2 | sort -k2,2) - | \
  cut -d' ' -f1 --complement | \
  parallel -j12 -C' ' --env sumstats '
    zgrep -w {4} ${sumstats}/EA/{1}.PHENO1.glm.linear.gz | \
    awk -vseqid={1} -vsnpid={2} -vsnpid38={3} -vrsid={4} -vprot={5} -vOFS="\t" "{print seqid,snpid,snpid38,rsid,prot,\$0}"
  '
) > ${INF}/ARIC/replication.tsv

(
  echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
  Rscript -e '
     INF <- Sys.getenv("INF")
     suppressMessages(library(dplyr))
     replication <- read.delim(file.path(INF,"ARIC","replication.tsv")) %>%
                    filter(P<=5e-8) %>%
                    mutate(prot_snpid=paste0(prot,"-",snpid),proxies="as sentinel",r2=1,
                           p=if_else(snpid!="chr2:102992675_C_T",as.character(P),"1.16196e-517")) %>%
                    left_join(within(read.delim(file.path(INF,"work","INF1.METAL")),{prot_snpid=paste0(prot,"-",MarkerName)})) %>%
                    left_join(select(pQTLtools::inf1,target.short,target,prot)) %>%
                    select(target.short,snpid,uniprot,rsid,cis.trans,proxies,r2,p,target)
     write.table(data.frame(replication,Source="Zhang et al. (2022)",PMID="",Comment=""),
                            col.names=FALSE,row.names=FALSE,sep="\t")
  '
) > ${INF}/ARIC/ARIC.tsv

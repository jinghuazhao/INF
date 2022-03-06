#!/usr/bin/bash

if [ ! -d ${INF}/ARIC ]; then mkdir ${INF}/ARIC; fi

export ARIC=~/rds/results/public/proteomics/ARIC

(
  cat <(echo -e "seqid\tsnpid\tsnpid38\trsid\tprot") | \
  paste - ${ARIC}/glm.linear.hdr
  cat <(echo pos38) <(cut -f2 ${INF}/work/INF1.b38) | \
  paste ${INF}/work/INF1.METAL - | \
  awk 'NR>1{snpid=$1;gsub(/chr[0-9]*:[0-9]*/,"",$1);print snpid,$0}' | \
  awk '{$2="chr"$5":"$23$2;print snpid,$1,$2,$3,$4,$21}' | \
  sort -k5,5 | \
  join -12 -25 <(sed '1d' ${ARIC}/seqid.txt | cut -f1,2 | sort -k2,2) - | \
  cut -d' ' -f1 --complement | \
  parallel -j12 -C' ' --env ARIC '
    zgrep -w {4} ${ARIC}/EA/{1}.PHENO1.glm.linear.gz | \
    awk -vseqid={1} -vsnpid={2} -vsnpid38={3} -vrsid={4} -vprot={5} -vOFS="\t" "{print seqid,snpid,snpid38,rsid,prot,\$0}"
  '
) > ${INF}/ARIC/replication.tsv

function region()
{
  (
    awk -v OFS='\t' '{print $0,"Prot","MarkerName","sentinel"}' ${ARIC}/glm.linear.hdr
    join -12 -21 <(sed '1d' ${ARIC}/seqid.txt | cut -f1-3 | sort -k2,2) \
                 <(Rscript -e '
                      INF <- Sys.getenv("INF")
                      suppressMessages(library(dplyr))
                      load(file.path(INF,"work","novel_data.rda"))
                      metal <- read.delim(file.path(INF,"work","INF1.METAL38")) %>%
                               select(prot,MarkerName,pos38)
                      novel_data <- left_join(novel_data,metal) %>%
                                    mutate(region=paste0(Chromosome,":",pos38-1e6,"-",pos38+1e6)) %>%
                                    select(uniprot,region,prot,MarkerName,rsid) %>%
                                    arrange(uniprot)
                      write.table(novel_data,row.names=FALSE,col.names=FALSE,quote=FALSE)
                   ') | \
    parallel -C' ' -j20 '
      tabix ${ARIC}/EA/{2}.PHENO1.glm.linear.gz {4} | awk -v prot={5} -vsnpid={6} -vrsid={7} -vOFS="\t" "{print \$0,prot,snpid,rsid}"
    '
  ) > ${INF}/ARIC/region38.tsv
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    rep38 <- read.delim(file.path(INF,"ARIC","region38.tsv")) %>%
             filter(P<=5e-8) %>%
             select(-A1_FREQ,-TEST,-OBS_CT,-T_STAT,-ERRCODE) %>%
             mutate(seqname=X.CHROM,start=as.integer(POS-1),end=POS) %>%
                    arrange(X.CHROM,POS)
    gr <- with(rep38,GenomicRanges::GRanges(seqnames=seqname,IRanges::IRanges(start,end,names=ID))) %>% unique()
    suppressMessages(library(rtracklayer))
    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- import.chain(path)
    seqlevelsStyle(gr) <- "UCSC"
    gr19 <- liftOver(gr,ch)
    df <- transmute(data.frame(gr),group_name=names(gr)) %>%
          left_join(data.frame(gr19)) %>%
          filter(!is.na(seqnames)) %>%
          mutate(pos37=end,Name=group_name) %>%
          select(Name,seqnames,pos37,group) %>%
          left_join(select(rep38,-X.CHROM,-seqname,-start,-end),by=c('Name'='ID')) %>%
          select(-group)
    write.table(subset(df,!is.na(Name)),file=file.path(INF,"ARIC","region.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
  '
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    region <- read.delim(file.path(INF,"ARIC","region.tsv")) %>%
              mutate(rsid=gsub("_[A-Z]*_[A-Z]*","",Name))
    key <- Sys.getenv("LDLINK_TOKEN")
    sentinels <- unique(region$sentinel)
    print(sentinels)
    blocks <- r <- list()
    check <- function(snps)
    {
      sentinel_and_snps <- c(sentinels[i],snps[grepl("^rs",snps)])
      r[[i]] <- LDlinkR::LDmatrix(snps=sentinel_and_snps,pop="EUR",r2d="r2",token=key)
      r2 <- subset(r[[i]],RS_number==sentinels[i])
      sel <- !is.na(r2) & r2>=0.8
      cat(sentinels[i],"\n")
      print(names(r2)[sel])
      print(r2[sel])
    }
    for(i in 1:length(sentinels))
    {
      blocks[[i]] <- subset(region,sentinel==sentinels[i])
      check(blocks[[i]]$rsid)
    }
  '
}

# > sentinels
# [1] "rs60094514" "rs7213460"
# rs60094514
# [1] "RS_number"  "rs60094514"
# [1] "rs60094514" "1"
#
# rs7213460
# [1] "RS_number" "rs7213460"
# [1] "rs7213460" "1"

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
                    select(target.short,rsid,uniprot,snpid,cis.trans,proxies,r2,p,target)
     write.table(data.frame(replication,Source="Zhang et al. (2022)",PMID="",Comment=""),
                            col.names=FALSE,row.names=FALSE,sep="\t")
  '
) > ${INF}/ARIC/ARIC.tsv

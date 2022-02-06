#!/usr/bin/bash

if [ ! -d ${INF}/AGES ]; then mkdir ${INF}/AGES; fi
cd ${INF}/AGES

export stables=https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx

function download()
{
R --no-save <<END
  INF <- Sys.getenv("INF")
  stables <- Sys.getenv("stables")
  st16 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, startRow=3)
  write.table(st16,file="AGES.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
  suppressMessages(library(dplyr))
  url1 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000"
  url2 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90087001-GCST90088000"
  url3 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90088001-GCST90089000"
  url4 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90089001-GCST90090000"
  url5 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90090001-GCST90091000"
  urls <- mutate(st16, tag=gsub("_",".",Study.tag),
                 trait=gsub("Serum levels of protein ","",Reported.trait),
                 Accession=Study.Accession,
                 SOMAMER_ID=paste0(trait,".",tag)) %>%
          left_join(pQTLtools::SomaLogic160410) %>%
          select(Accession,Summary.statistics.file,SOMAMER_ID,UniProt,Target) %>%
          filter(UniProt %in% pQTLtools::inf1$uniprot) %>%
          distinct() %>%
          mutate(url=case_when(
                                 Accession >= "GCST90086001" & Accession <= "GCST90087000" ~ url1,
                                 Accession >= "GCST90087001" & Accession <= "GCST90088000" ~ url2,
                                 Accession >= "GCST90088001" & Accession <= "GCST90089000" ~ url3,
                                 Accession >= "GCST90089001" & Accession <= "GCST90090000" ~ url4,
                                 Accession >= "GCST90090001" & Accession <= "GCST90091000" ~ url5,
                                 TRUE ~ as.character(Accession)
                              ),
                 src=file.path(url,Accession),cmd=paste("lftp -c mirror",src))
  write.table(select(urls,Accession,SOMAMER_ID,UniProt,Target) %>% left_join(select(pQTLtools::inf1,-target),by=c('UniProt'='uniprot')),
              file=file.path(INF,"AGES","links.txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(select(urls,cmd),file=file.path(INF,"AGES","lftp.sh"),quote=FALSE,col.names=FALSE,row.names=FALSE)
END

cd ~/rds/results/public/proteomics/AGES
bash ~/INF/AGES/lftp.sh
ls GC*/*tsv | parallel -C' ' -j15 'bgzip {}'
ls GC*/*gz | parallel -C' ' -j15 'tabix -S1 -s3 -b4 -e4 -f {}'
cd -

(
  gunzip -c ~/rds/results/public/proteomics/AGES/GCST*/*tsv.gz | \
  head -1 | \
  awk -vOFS="\t" '{print "Protein","Comment",$0}'
  join -13 -22 <(cut -f1,2,3 ${INF}/work/INF1.METAL | sed '1d' | sort -k3,3) \
               <(sed '1d' ${INF}/AGES/links.txt | grep -v BDNF | cut -f1,5,7 | sort -k2,2 ) | \
  parallel -j15 -C' ' '
    zgrep -w {3} ~/rds/results/public/proteomics/AGES/{4}/{4}_buildGRCh37.tsv.gz | \
    awk -vOFS="\t" -vprot="{1}\t{4}" "{print prot,\$0}"
  '
) > ${INF}/AGES/INF1.replication.txt

cat <(awk 'NR==1' ${INF}/AGES/INF1.replication.txt) \
    <(awk 'NR==1||$4<=5e-8' ${INF}/AGES/INF1.replication.txt | \
      sort -k1,1 -k2,2 -k4,4g | \
      awk 'a[$1$3]++==0') > ${INF}/AGES/INF1.replication

Rscript -e '
  INF <- Sys.getenv("INF")
  suppressMessages(library(dplyr))
  replication <- read.delim(file.path(INF,"AGES","INF1.replication")) %>%
                 select(Protein,variant_id,p_value,Comment) %>%
                 rename(p=p_value)
  metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
           select(MarkerName,rsid,prot,uniprot,cis.trans) %>%
           left_join(replication,by=c('prot'='Protein','rsid'='variant_id')) %>%
           filter(!is.na(p)) %>%
           left_join(pQTLtools::inf1[c("prot","target.short","target")]) %>%
           rename(Protein=target.short,Sentinels=rsid,UniProt=uniprot,SNPid=MarkerName,"cis/trans"=cis.trans,TargetFullName=target) %>%
           mutate(Source="Gudjonsson, et al. (2022)",PMID=NA,Proxies="as sentinels",r2=1,p=as.character(p)) %>%
           select(Protein,Sentinels,UniProt,SNPid,"cis/trans",Proxies,r2,p,TargetFullName,Source,PMID,Comment)
  write.table(metal,file=file.path(INF,"AGES","AGES.tsv"),row.names=FALSE,sep="\t")
'

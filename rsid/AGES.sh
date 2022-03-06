#!/usr/bin/bash

export AGES=~/rds/results/public/proteomics/AGES

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

cd ${AGES}
bash ~/INF/AGES/lftp.sh
ls GC*/*tsv | parallel -C' ' -j15 'bgzip {}'
ls GC*/*gz | parallel -C' ' -j15 'tabix -S1 -s3 -b4 -e4 -f {}'
cd -

(
  awk -vOFS="\t" '{print "Protein","Comment",$0}' ${AGES}/AGES.hdr
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

function region()
{
  (
    awk -vOFS="\t" '{print $0,"Protein","MarkerName","sentinel"}' ${AGES}/AGES.hdr
    join -11 -22 <(Rscript -e '
                     INF <- Sys.getenv("INF")
                     suppressMessages(library(dplyr))
                     load(file.path(INF,"work","novel_data.rda"))
                     novel_data <- novel_data %>%
                                   mutate(region=paste0(Chromosome,":",Position-1e6,"-",Position+1e6)) %>%
                                          select(prot,MarkerName,rsid,region)
                     write.table(novel_data,row.names=FALSE,col.names=FALSE,quote=FALSE)
                   ' | sort -k1,1) \
                 <(sed '1d' ${INF}/AGES/links.txt | grep -v BDNF | cut -f1,5,7 | sort -k2,2 ) | \
    parallel -j15 -C' ' '
      tabix ~/rds/results/public/proteomics/AGES/{5}/{5}_buildGRCh37.tsv.gz {4} | \
      awk -vOFS="\t" -vprot={1} -vsnpid={2} -vrsid={3} "{print \$0,prot,snpid,rsid}"
    '
  ) > ${INF}/AGES/region.tsv
  Rscript -e '
    INF <- Sys.getenv("INF")
    suppressMessages(library(dplyr))
    region <- read.delim(file.path(INF,"AGES","region.tsv")) %>%
              mutate(rsid=gsub("_[A-Z]*_[A-Z]*","",variant_id))
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

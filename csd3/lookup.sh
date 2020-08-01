#!/usr/bin/bash

function check()
# perform check over an individual variant
{
  grep -w $1 work/INF1.merge.rsid
  grep -w $1 work/INF1.merge.cis.vs.trans
  grep -w $1 work/INF1.merge.cis.vs.trans | cut -d' ' -f2 | sed 's/"//g'
  export prot=$(grep -w $1 work/INF1.merge.cis.vs.trans | cut -d' ' -f2)
  grep ${prot} doc/olink.inf.panel.annot.tsv
}

function Sun()
{
  grep -w $1 work/ps/pQTL.Sun-B_pQTL_EUR_2017
  check $1
}

function Folkersen()
{
  grep -w $1 work/ps/pQTL.Folkersen-L_Proteins_EUR_2017
  check $1
}

function Suhre()
{
  grep -w $1 work/ps/pQTL.Suhre-K_pQTL_EUR_2017
  check $1
}

function pQTL()
{
  grep -w $1 work/ps/pQTL.pQTL_2017
  check $1
}

for rsid in $(sed '1d' work/ps/pQTL.Sun-B_pQTL_EUR_2017 | awk '{print $3}' | uniq)
do
  echo --- ${rsid} ---
  Sun ${rsid}
  echo
done

function lookup()
# extended lookup for Olink-SomaLogic signals
{
# Olink sentinel thread over SomaLogic
  join <(awk 'NR>1{print $5,"chr" $8 ":" $9}' work/INF1.merge | sort -k1,1 -k2,2) <(sort -k1,1 work/inf1.tmp) | \
  parallel -C' ' '
    echo {1}+{2}+{3}
    grep {2} pQTL.Sun-B_pQTL_EUR_2017 | \
    grep {3}
  ' > pQTL.Sun.log
# Olink + SomaLogic signal overlap
  export s=work/pQTL_2018.txt.gz
  awk -v OFS='\t' '{gsub(/"/,"");print $1,$2,$3":"$4,$5}' work/INF1.merge.cis.vs.trans | \
  sort -k1,1 | \
  join -25 -t$'\t' - <(sed '1d' INTERVAL_box.tsv | grep -v BDNF | sort -t$'\t' -k5,5) | \
  cut -f1,2,3,4,10 | \
  sort -t$'\t' -k5,5 | \
  join -t$'\t' -j5 - <(zcat ${s} | sed '1d' | awk -v FS='\t' '/29875488/ && $12 <= 1.5e-11' | sort -t $'\t' -k5,5) > Olink+SomaLogic.ps
  R --no-save -q <<\ \ END
    options(width=160)
    library(dplyr)
    # by SNP
    os <- read.delim("Olink+SomaLogic.ps", header=FALSE, col.names=c("trait","uniprot","prot","chrpos","SNP","rsid","hg19_coordinates",
                     "a1","a2","efo","study","pmid","ancestry","beta","se","p","direction","n","n_studies","unit","dataset"))
    v <- c("trait","uniprot","prot","chrpos","SNP","rsid","hg19_coordinates","a1","a2","beta","se","p")
    m <- os[v] %>% group_by(prot,rsid) %>% slice(which.min(p))
    # by gene
    library(gap)
    g <- with(subset(inf1,gene!="BDNF"),gene)
    batches <- split(g, ceiling(seq_along(g)/10))
    s <- t <- list()
    for(i in 1:length(batches))
    {
      cat("Block ",i,batches[[i]],"\n")
      q <- phenoscanner(genequery=batches[[i]], catalogue="pQTL", proxies="EUR", pvalue=1.5e-11, r2=0.7, build=37)
      s[[i]] <- with(q,genes)
      t[[i]] <- with(q,results)
    }
    r <- list(genes=do.call(rbind,s),results=within(do.call(rbind,t),{
       ref_a1 <- a1
       ref_a2 <- a2
       swap <- a1 > a2
       a1[swap] <- ref_a2[swap]
       a2[swap] <- ref_a1[swap]
       snpid <- paste0(hg19_coordinates,"_",a1,"_",a2)
       p <- as.numeric(p)
    }))
    setdiff(g,with(r,genes)[["gene"]])
    # Sun, et al. (2018)
    vars <- c("gene", "rsid", "snpid", "hg19_coordinates", "a1", "a2", "ensembl", "hgnc", "pmid", "beta", "se", "p")
    sun <- subset(with(r,results[vars]),pmid==29875488&gene==hgnc)
    length(table(sun$gene))
    length(table(sun$snpid))
    m <- sun %>% group_by(gene,rsid) %>% slice(which.min(p))
    for(g in unique(with(m,gene))) print(subset(m,gene==g))
    nosig <- read.table("work/INF1.merge.nosig",as.is=TRUE,col.names=c("prot","uniprot"))
    nosig <- subset(inf1,uniprot%in%nosig[["uniprot"]])
    SomaLogic_yes_olink_no <- subset(sun,gene%in%nosig[["gene"]])
  END
}

Sun

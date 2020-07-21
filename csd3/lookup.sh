#!/usr/bin/bash

function check()
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

function Sun()
{
  join <(awk 'NR>1{print $5,"chr" $8 ":" $9}' work/INF1.merge | sort -k1,1 -k2,2) <(sort -k1,1 work/inf1.tmp) | \
  parallel -C' ' '
    echo {1}+{2}+{3}
    grep {2} pQTL.Sun-B_pQTL_EUR_2017 | \
    grep {3}
  ' > pQTL.Sun.log
# Olink + SomaLogic
  cut -f2 work/INF1.merge.nosig | \
  grep -f - -v work/inf1.tmp | \
  sort -k2,2 | \
  join -12 -25 -t$'\t' - <(sed '1d' INTERVAL_box.tsv | sort -t$'\t' -k5,5) | \
  cut -f1,2,8 | \
  sort -t$'\t' -k3,3 | \
  join -t$'\t' -13 -25 - <(zcat work/pQTL_2018.txt.gz | sed '1d' | awk -v FS='\t' '/29875488/ && $12 <= 1.5e-11' | sort -t $'\t' -k5,5) > Olink+SomaLogic.ps
# 
  awk -v FS='\t' '!/BDNF/ && NR > 1 {
    gsub(/"/,"",$7);if ($3=="Q8NF90") $7="FGF5"; else if ($3=="Q8WWJ7") $7="CD6"; print $7
  }' doc/olink.inf.panel.annot.tsv > inf.genes
  R --no-save -q <<\ \ END
    options(width=160)
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
    library(dplyr)
    m <- sun %>% group_by(gene,rsid) %>% slice(which.min(p))
    for(g in unique(with(m,gene))) print(subset(m,gene==g))
    nosig <- read.table("work/INF1.merge.nosig",as.is=TRUE,col.names=c("prot","uniprot"))
    nosig <- subset(inf1,uniprot%in%nosig[["uniprot"]])
    SomaLogic_yes_olink_no <- subset(sun,gene%in%nosig[["gene"]])
    xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
    t4 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:31), rows=c(5:1986))
    t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:19), rows=c(3:2746))
    t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:20), rows=c(3:167))
  END
}

Sun

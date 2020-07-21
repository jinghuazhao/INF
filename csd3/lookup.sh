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
    genes <- scan("inf.genes","")
    library(phenoscanner)
    ps <- phenoscanner(genequery=genes[1:10],catalogue="pQTL",pvalue=1.5e-11)
    vars <- c("gene", "rsid", "hg19_coordinates", "a1", "a2", "ensembl", "hgnc", "pmid", "beta", "se", "p")
    bys <- c("gene","rsid","p")
    r <- subset(with(ps,results[vars]),pmid==29875488&gene==hgnc)
    attach(r)
    m <- aggregate(r[bys],by=list(gene,rsid),FUN=min,na.rm=TRUE)
    detach(r)
    merge(r,m[bys],by=c("gene","rsid","p"))
    # Sun, et al. (2018)
    xlsx <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"
    t4 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:31), rows=c(5:1986))
    t5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:19), rows=c(3:2746))
    t6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:20), rows=c(3:167))
  END
}

Sun

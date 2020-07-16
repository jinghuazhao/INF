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
#   grep -H -w {3} INTERVAL_box.tsv
  ' > pQTL.Sun.log
  join -12 -25 -t$'\t' <(sort -k2,2 work/inf1.tmp) <(sed '1d' INTERVAL_box.tsv | sort -k5,5) > Olink+SomaLogic.list
  cut -f1 Olink+SomaLogic.list | \
  grep -v P23560 | \
  uniq | \
  wc -l
  join -17 -25 -t$'\t' <(sed '1d' INTERVAL_box.tsv | sort -k7,7) <(gunzip -c work/pQTL_2018.txt.gz | sed '1d' | sort -k5,5) | \
  sort -k6,6 | \
  join -16 -t$'\t' - <(cut -f1 Olink+SomaLogic.list) > Olink+SomaLogic.ps
}

Sun

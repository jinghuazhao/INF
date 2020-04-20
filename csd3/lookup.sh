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
  join <(awk 'NR>1{print $5,"chr" $8 ":" $9}' work/INF1.merge | sort -k1,1 -k2,2) work/inf1.tmp | \
  parallel -C' ' '
    echo {1} {2} {3}
    grep -w {3} pQTL.Sun-B_pQTL_EUR_2017 | grep {2}
    grep -H -w {3} INTERVAL_box.tsv
  '
}

Sun

function Olink()
{
  export OLINK=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output
  ls $OLINK > olink.list
  join -11 -22 \
       <(awk 'NR>1{print $5,$8 ":" $9}' work/INF1.merge | sort -k1,1) \
       <(ls $OLINK/*gz  | sed 's/___/ /g;s/_chr_merged.gz\*//g;s///g;s///g;s///g' | sort -k2,2) | \
  awk '{
     gsub(/chr/,"",$2);
     split($2,a,":");
     chr=a[1];
     pos=a[2];
     print $1,pos,$3,chr
  }' | \
  parallel -C' ' '
    echo {1} {2}
    zgrep -H -w {2} {3}___{1}_chr_merged.gz | \
    awk -vchr={4} "(\$3==chr)"
  '
}

(
  gunzip -c /rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output/INTERVAL_cvd3_SELP___P16109_chr_merged.gz | \
  awk 'NR==1{print "UniProt","pos",$2,$22,$24,$25}'
  Olink | \
  awk '{if(NF==2) printf $1,$2," "; else print $2,$22,$24,$25}'
) > olink.overlap


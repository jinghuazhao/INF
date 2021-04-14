#!/usr/bin/bash

if [ ! -d ${INF}/work/HGI ]; then mkdir ${INF}/work/HGI; fi

function _outcome_r4()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
  grep -e chr3 -e chr9 work/INF1.merge-rsid | \
  cut -f5,6,8,9 --output-delimiter=' ' | \
  parallel -C' ' '
     gunzip -c ${HGI}/COVID19_HGI_C2_ALL_leave_23andme_20201020.b37.txt.gz | \
     awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > work/HGI/{1}-{2}
     cat ${HGI}/COVID19_HGI_C2_ALL_20201020.b37_1.0E-5.txt | \
     awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > work/HGI/{1}-{2}_1e-5
  '
}

function _outcome_r5()
{
  export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107
  for trait in A2 B2 C2
  do
    grep -e chr3 -e chr9 ${INF}/work/INF1.merge-rsid | \
    cut -f5,6,8,9 --output-delimiter=' ' | \
    parallel -C' ' '
       gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
       awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > ${INF}/HGI/${trait}-{1}-{2}
    '
  done
}
_outcome_r5

function _exposure()
{
  (
    echo SNP Phenotype Allele1 Allele2 EAF Effect StdErr P N
    grep -e chr3 -e chr9 ${INF}/work/INF1.METAL | \
    cut -f2,3,6,7,8-11,17,21 | \
    awk '{print $1,$2,toupper($3),toupper($4),$5,$6,$7,10^$8,$9}'
  ) > ${INF}/HGI/INF.ins
  for trait in A2 B2 C2
  do
    grep -e chr3 -e chr9 ${INF}/work/INF1.merge-rsid | \
    cut -f5,6,8,9 --output-delimiter=' ' | \
    parallel -j5 -C' ' '
      gunzip -c ${INF}/METAL/{1}-1.tbl.gz | \
      cut -f1-6,10-12,18 | \
      awk -vchr={3} -vpos={4} -vM=1e6 -vlogp=-5.45131 "\$1==chr && \$2>=pos-M && \$2 <= pos+M && \$9<=logp" > work/HGI/${trait}-{1}-{2}.mri
      (
        echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
        awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" ${INF}/HGI/{1}-{2}.mri | \
        sort -k3,3 | \
        join -23 ${INF}/work/INTERVAL.rsid - | \
        awk -v prot={1} "{\$1=prot;print}" | \
        tr " " "\t"
      ) | gzip -f > ${INF}/HGI/-${trait}{1}-{2}.mrx
    '
  done
}
_exposure

module load gcc/6
R --no-save < rsid/HGI.R

# grep -e chr3 -e chr9 work/INF1.merge-rsid | cut -f6 | zgrep -f - -w $C2
ls *-rs* | sed 's/-/ /g' | \
parallel -C' ' '
    echo {1} - {2};
    grep -w {2} {1}-{2};
    grep -w -e {1} work/INF1.merge.cis.vs.trans-rsid | grep -w {2};
    grep -w -e {1} INF.ins | grep -e {2}
'

#!/usr/bin/bash

if [ ! -d ${INF}/HGI/mr ]; then mkdir ${INF}/HGI/mr; fi

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
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    parallel --env trait -C' ' '
       gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
       awk -v M=1e6 -v chr={3} -v pos={4} "NR==1||(\$1==chr && \$2>=pos-M && \$2<=pos+M)" > ${INF}/HGI/mr/${trait}-{1}-{2}
    '
  done
}
_outcome_r5

function _ins()
(
  echo SNP Phenotype Allele1 Allele2 EAF Effect StdErr P N
  sed '1d' ${INF}/work/INF1.METAL | \
  cut -f2,3,6,7,8-11,17,21 | \
  awk '{print $1,$2,toupper($3),toupper($4),$5,$6,$7,10^$8,$9}'
) > ${INF}/HGI/mr/INF.ins
_ins

function _exposure_vcf()
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    awk -vM=1e6 '{start=$4-M;if(start<0) start=1;end=$4+M;print $1,$2,$3":"start"-"end}' | \
    parallel -j15 --env INF --env trait -C' ' '
    (
      echo -e "prot\trsid\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN\tChromosome\tPosition"
      bcftools query -f "{1}\t%ID\t%ALT\t%REF\t%AF\t[%ES]\t[%SE]\t[%LP]\t[%SS]\t%CHROM\t%POS\n" \
                     -r {3} ${INF}/METAL/gwas2vcf/{1}.vcf.gz | \
      awk "a[\$2]++==0"
    ) | \
    gzip -f > ${INF}/HGI/mr/${trait}-{1}-{2}.gz
    '
  done
_exposure_vcf

module load gcc/6

for trait in A2 B2 C2
do
  export trait=${trait}
  sed '1d' ${INF}/HGI/mr/INF.ins | \
  parallel -C' ' '
    export prot={2}; export rsid={1};
    export Allele1={3}; export Allele2={4}; export EAF={5}
    export Effect={6}; export StdErr={7}; export P={8}; export N={9}
    R --no-save < ${INF}/rsid/HGI.R
  '
done

function _exposure_tsv()
{
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    awk -vM=1e6 '{start=$4-M;if(start<0) start=1;end=$4+M;print $1,$2,$3":"start"-"end}' | \
    parallel -j15 --env INF --env trait -C' ' '
    (
      echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP\tN"
      tabix ${INF}/METAL/gwas2vcf/{1}.tsv.gz {3} | \
      awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" | \
      sort -k3,3 | \
      join -23 ${INF}/work/INTERVAL.rsid - | \
      awk -v prot={1} "a[\$2]++==0{\$1=prot;print}" | \
      tr " " "\t"
    ) | gzip -f > ${INF}/HGI/mr/${trait}-{1}-{2}.tsv.gz
    '
  done
}
_exposure_tsv

function _exposure_METAL()
{
  for trait in A2 B2 C2
  do
    export trait=${trait}
    cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge-rsid | \
    sed '1d' | \
    parallel -j15 --env trait -C' ' '
      gunzip -c ${INF}/METAL/{1}-1.tbl.gz | \
      cut -f1-6,10-12,18 | \
      awk -vchr={3} -vpos={4} -vM=1e6 "\$1==chr && \$2>=pos-M && \$2 <= pos+M" > ${INF}/HGI/mr/${trait}-{1}-{2}.mri
      (
        echo -e "prot\trsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
        awk "{\$4=toupper(\$4);\$5=toupper(\$5);print}" ${INF}/HGI/mr/${trait}-{1}-{2}.mri | \
        sort -k3,3 | \
        join -23 ${INF}/work/INTERVAL.rsid - | \
        awk -v prot={1} "a[\$2]++==0{\$1=prot;print}" | \
        tr " " "\t"
      ) | gzip -f > ${INF}/HGI/mr/${trait}-{1}-{2}.mrx
    '
  done
}
_exposure_METAL

function MR_collect()
{
  awk -vOFS="\t" '/Inverse/ && $NF<0.05/180 {print FILENAME,$9,$10,$11}' ${INF}/HGI/mr/MR*result*
}

function pqtlMR_collect()
{
  awk -vOFS="\t" '$NF<0.05/180 {print FILENAME,$8,$9,$10}' pqtl*result* | \
  awk '{gsub(/pqtlMR-|-result.txt/,"",$1)};1' | \
  awk -vOFS="\t" '{gsub(/-/,"\t",$1)};1' | \
  xsel -i
}

function etc()
{
# grep -e chr3 -e chr9 work/INF1.merge-rsid | cut -f6 | zgrep -f - -w $C2
  ls *-rs* | sed 's/-/ /g' | \
  parallel -C' ' '
      echo {1} - {2};
      grep -w {2} {1}-{2};
      grep -w -e {1} work/INF1.merge.cis.vs.trans-rsid | grep -w {2};
      grep -w -e {1} INF.ins | grep -e {2}
  '
  ls ${INF}/HGI/mr/no-clumping/MR*result* | \
  sed 's/-/ /g;' | \
  cut -d' ' -f5-7 | \
  parallel -C' ' 'echo {1} {2} {3}'
}

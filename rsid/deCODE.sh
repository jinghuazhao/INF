#!/usr/bin/bash

export proteomics_results=~/rds/results/public/proteomics/deCODE/full_downloads
export v4=SomaLogicv4.tsv

if [ ! -d ${INF}/deCODE ]; then mkdir -p ${INF}/deCODE; fi

function replication()
(
  gunzip -c ${proteomics_results}/10000_28_CRYBB2_CRBB2.txt.gz | \
  head -1
  join -j2 <(cut -f1,3,6 --complement --output-delimiter=' ' ${INF}/Fenland/${v4} | sed 's/-/_/' | sort -k2,2) \
           <(cut -f2,4,5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | awk '{print $2":"$3,$4,$1}' | sort -k2,2) | \
  parallel -C' ' -j20 --env proteomics_results '
    tabix ${all} {4} | zgrep -w {2} | grep -w {5}
  '
) > ${INF}/deCODE/replication.tsv


function olink_inf_idx()
{
  module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7
  cd ${INF}/deCODE
  ls ${proteomics_results}/full_downloads/*md5sum | xargs -l basename -s .md5sum | \
  parallel -C' ' '
    gunzip -c 
    tabix -f -S1 -s1 -b2 -e2 {}.gz
  '
  cd -
}

function deCODE()
(
echo -e "Protein\tSentinels\tUniProt\tSNPid\tcis/trans\tProxies\tr2\tp\tTarget Full Name\tSource\tPMID\tComment"
join -12 -21 <(awk 'NR>1 && $14<=5e-8' ${INF}/deCODE/replication.tsv | cut -f4,11,14 | sort -k2,2) \
             <(cut -f1,3,5,6 --complement --output-delimiter=' ' ${INF}/Fenland/${v4} | sed 's/-/_/' | sort -k1,1) | \
awk '{print $4"-"$2,$0}' | \
sort -k1,1 | \
join - <(awk 'NR>1 {print $20"-"$1,$3,$2,$21}' ${INF}/work/INF1.METAL | sort -k1,1) | \
awk 'a[$1]++==0' | \
cut -d' ' -f1,2 --complement | \
sort -k4,4 | \
tr ' ' '\t'| \
join -24 <(Rscript -e 'write.table(pQTLtools::inf1[c("prot","target.short")],col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
           sort -t$'\t' -k1,1) - | \
awk -vOFS='\t' '{print $2,$6,$5,$3,$7,"as sentinels",1,$4}' | \
sort -t$'\t' -k1,1 | \
join -t$'\t' - \
             <(Rscript -e 'write.table(data.frame(pQTLtools::inf1[c("target.short","target")],Source="Ferkingstad et al. (2021)",PMID="",Comment=""),
                                       col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")' | \
               sort -t$'\t' -k1,1)
)

deCODE > ${INF}/deCODE/deCODE.tsv

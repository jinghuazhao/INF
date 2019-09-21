# 21-9-2019 JHZ

export TMPDIR=/rds/user/jhz22/hpc-work/work
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export tag=_nold

for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do

echo $p
awk -v OFS="\t" '
{
   if (NR==1) print "MarkerName","EPACTS","CHR","POS","PVALUE";
   else {
     gsub(/chr/,"",$1)
     EPACTS=$1 ":" $3 "_" toupper($5) "/" toupper($6)
     print $4,EPACTS,$1,$3,exp($13)
   }
}' ${INF}/sentinels/${p}${tag}.p > work/${p}.p

swiss --assoc ${INF}/work/${p}.p \
      --variant-col EPACTS \
      --chrom-col CHR \
      --pos-col POS \
      --pval-col PVALUE \
      --include-cols MarkerName --trait ${p} \
      --dist-clump --clump-dist 1000000 --clump-p 5e-10 --out ${INF}/work/${p}
done

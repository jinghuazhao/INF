export swiss_data=$HOME/.local/share/swiss/data

awk -v OFS="\t" '
{
   if (NR==1) print "MarkerName","CHR","POS","PVALUE";
   else {
     gsub(/chr/,"",$1)
     print $4,$1,$3,exp($13)
   }
}' sentinels/IL.6_nold.p > work/IL.6.p

swiss --assoc work/IL.6.p \
      --build hg19 \
      --ld-clump-source 1000G_2014-11_EUR \
      --variant-col MarkerName \
      --chrom-col CHR \
      --pos-col POS \
      --pval-col PVALUE \
      --dist-clump --clump-dist 1000000 --clump-p 5e-10 --out test

#     --ld-gwas-source 1000G_2014-11_EUR \
#     --gwas-cat data/gwascat_ebi_GRCh37p13.tab \



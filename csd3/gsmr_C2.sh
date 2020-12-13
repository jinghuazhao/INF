export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
(
  echo "SNP A1 A2 freq b se p N"
  awk '
  {
    CHR=$1
    POS=$2
    a1=$4
    a2=$3
    if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
    else snpid="chr" CHR ":" POS "_" a1 "_" a2
    if (NR>1) print snpid, a1, a2, $12, $7, $8, $9, $11
  }' ${HGI}/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37_1.0E-5.txt
) | \
gzip -f > work/HGI/gsmr_C2.txt.gz

(
  echo "SNP A1 A2 freq b se p N"
  gunzip -c $HGI/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37.txt.gz | \
  awk '
  {
    CHR=$1
    POS=$2
    a1=$4
    a2=$3
    if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
    else snpid="chr" CHR ":" POS "_" a1 "_" a2
    if (NR>1) print snpid, a1, a2, $12, $7, $8, $9, $11
  }'
) | \
awk 'a[$1]++==0' | \
gzip -f > work/HGI/gsmr_C2.txt.gz

#  head -1 ${HGI}/COVID19_HGI_C2_ALL_20201020.b37_1.0E-5.txt | sed 's/\t/\n/g' | awk '{print "#" NR, $1}'
#1 #CHR
#2 POS
#3 REF
#4 ALT
#5 SNP
#6 all_meta_N
#7 all_inv_var_meta_beta
#8 all_inv_var_meta_sebeta
#9 all_inv_var_meta_p
#10 all_inv_var_het_p
#11 all_meta_sample_N
#12 all_meta_AF
#13 rsid

if [ -f work/HGI/INF1_C2.gsmr ]; then rm work/HGI/INF1_C2.gsmr; fi
(
  cat work/HGI/gsmr_C2*.gsmr | \
  head -1
  ls work/HGI/gsmr_C2*gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > work/HGI/INF1_C2.gsmr

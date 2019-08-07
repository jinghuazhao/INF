# 7-8-2019 JHZ

(
  echo "SNP A1 A2 freq b se p N"
  unzip -p -c cad.additive.Oct2015.pub.zip | \
  awk '
  {
    CHR=$2
    POS=$3
    a1=$4
    a2=$5
    if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
    else snpid="chr" CHR ":" POS "_" a1 "_" a2
    if (substr($1,1,2)=="rs"||substr($1,1,3)=="chr"||substr($1,1,3)=="MER")
       print snpid, a1, a2, $6, $9, $10, $11, 185000
  }'
) > gsmr_CAD.txt

# unzip -p -c CAD/cad.additive.Oct2015.pub.zip | head -1 | sed 's/\t/\n/g' | awk '{print "#" NR, $1}'
#1 markername
#2 chr
#3 bp_hg19
#4 effect_allele
#5 noneffect_allele
#6 effect_allele_freq
#7 median_info
#8 model
#9 beta
#10 se_dgc
#11 p_dgc
#12 het_pvalue
#13 n_studies


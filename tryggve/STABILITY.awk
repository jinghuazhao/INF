# 19-3-2019 JHZ

{
  if (NR==1||($9>=1e-6 && $14>0.4)) print $1,$2,$3,$4,N,$6,$7,$8,$10,$11,$12,$13,$14,$15
}

# gunzip -c  work/STABILITY.4E.BP1.gz | head -1 | sed 's/ /\n/g' | awk '{print "#" NR, $1}'

#1 SNPID
#2 CHR
#3 POS
#4 STRAND
#5 N
#6 EFFECT_ALLELE
#7 REFERENCE_ALLELE
#8 CODE_ALL_FQ
#9 HWE_PVAL
#10 BETA
#11 SE
#12 PVAL
#13 RSQ
#14 RSQ_IMP
#15 IMP

{
   OFS="\t"
   if (NR>1)
   {
     CHR=$2
     POS=$3
     a1=$6
     a2=$7
     if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
     else snpid="chr" CHR ":" POS "_" a1 "_" a2
     $1=snpid
   }
   print
}

# 1 SNPID
# 2 CHR
# 3 POS
# 4 STRAND
# 5 N
# 6 EFFECT_ALLELE
# 7 REFERENCE_ALLELE
# 8 CODE_ALL_FQ
# 9 BETA
# 10 SE
# 11 PVAL
# 12 RSQ
# 13 RSQ_IMP
# 14 IMP

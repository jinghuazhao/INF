# 1-11-2018 JHZ

{
  OFS="\t"
  if (NR==1) print  "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  else {
    CHR=$2
    POS=$1
    if(substr($3,1,2)=="rs") SNPID=$3; else SNPID="chr" CHR ":" POS
    STRAND="NA"
    N=$8
    EFFECT_ALLELE=$7
    REFERENCE_ALLELE=$6
    CODE_ALL_FQ=$5
    BETA=$9
    SE=$10
    PVAL=$12
    RSQ="NA"
    RSQ_IMP="NA"
    IMP="NA"
    print  SNPID, CHR, POS, STRAND, N, EFFECT_ALLELE, REFERENCE_ALLELE, CODE_ALL_FQ, BETA, SE, PVAL, RSQ, RSQ_IMP, IMP
  }
}

#1 BP
#2 CHR
#3 SNP
#4 HWE
#5 MAF
#6 A1
#7 A2
#8 N
#9 BETA
#10 SE
#11 CHI2
#12 P
#13 P_GC

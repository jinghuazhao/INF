{
  OFS="\t"
  if (NR==1) print  "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  else {
    SNPID=$2
    CHR=$3
    sub(/^0/,"",CHR)
    POS=$4
    STRAND="+"
    N_AA=$14
    N_AB=$15
    N_BB=$16
    N=N_AA+N_AB+N_BB
    EFFECT_ALLELE=$6
    REFERENCE_ALLELE=$5
    CODE_ALL_FQ=(N_AA+N_AB*0.5)/N
    BETA=$24
    SE=$25
    PVAL=$22
    RSQ=$9
    RSQ_IMP=$23
    IMP=0
    print  SNPID, CHR, POS, STRAND, N, EFFECT_ALLELE, REFERENCE_ALLELE, CODE_ALL_FQ, BETA, SE, PVAL, RSQ, RSQ_IMP, IMP
  }
}

# to check, strand
# SNPTEST output
# h <- read.table("h", as.is=TRUE, header=TRUE)

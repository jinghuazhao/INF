# 19-2-2019 JHZ

{
  OFS="\t"
  if (NR==1) print  "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  else {
    CHR=$1
    POS=$2
    SNPID=$17
    STRAND="NA"
    N=$11
    EFFECT_ALLELE=$6
    REFERENCE_ALLELE=$4
    CODE_ALL_FQ=$9
    BETA=$12
    SE=$13
    PVAL=$15
    RSQ="NA"
    RSQ_IMP="NA"
    IMP="NA"
    print  SNPID, CHR, POS, STRAND, N, EFFECT_ALLELE, REFERENCE_ALLELE, CODE_ALL_FQ, BETA, SE, PVAL, RSQ, RSQ_IMP, IMP
  }
}

#1 #CHROM
#2 POS
#3 ID
#4 REF
#5 ALT
#6 A1
#7 A1_CT
#8 ALLELE_CT
#9 A1_FREQ
#10 TEST
#11 OBS_CT
#12 BETA
#13 SE
#14 T_STAT
#15 P
#16 ID:REF_ALT
#17 new_SNP

# head -1 /data/andmala/biofinder_inf/rsannot_runGwas_plasmaImp.FGF23_zre_INFI.glm.linear | \
# awk '{gsub(/\t/, "\n",$0)};1'| \
# awk '{print "#" NR, $1}'

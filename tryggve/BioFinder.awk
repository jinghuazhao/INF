# 19-2-2019 JHZ

{
  OFS="\t"
  if (NR==1) print  "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  else {
    split($1,a,":")
    CHR=a[1]
    POS=$4
    if (substr($2,1,2)=="rs") SNPID=$2; else SNPID="chr" CHR ":" POS
    STRAND="NA"
    N_AA=$14
    N_AB=$15
    N_BB=$16
    N=N_AA+N_AB+N_BB
    EFFECT_ALLELE=$6
    REFERENCE_ALLELE=$5
    CODE_ALL_FQ=(N_BB+N_AB*0.5)/N
    BETA=$24
    SE=$25
    PVAL=$22
    RSQ="NA"
    RSQ_IMP=$23
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

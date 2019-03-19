# 4-3-2019 JHZ

{
  OFS="\t"
  if (NR>1) {
     OFS="\t"
     CHR=$3
     sub(/^0/,"",CHR)
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
     RSQ_IMP=$9
     IMP="NA"
     a1=EFFECT_ALLELE
     a2=REFERENCE_ALLELE
     if (a1>a2) SNPID="chr" CHR ":" POS "_" a2 "_" a1;
     else SNPID="chr" CHR ":" POS "_" a1 "_" a2
     print  SNPID, CHR, POS, STRAND, N, EFFECT_ALLELE, REFERENCE_ALLELE, CODE_ALL_FQ, BETA, SE, PVAL, RSQ, RSQ_IMP, IMP
  }
}

#1 alternate_ids
#2 rsid, which treats as missing (.) for non-rsid
#3 chromosome
#4 position
#5 alleleA
#6 alleleB
#7 index
#8 average_maximum_posterior_call
#9 info
#10 cohort_1_AA
#11 cohort_1_AB
#12 cohort_1_BB
#13 cohort_1_NULL
#14 all_AA
#15 all_AB
#16 all_BB
#17 all_NULL
#18 all_total
#19 all_maf
#20 missing_data_proportion
#21 cohort_1_hwe
#22 frequentist_add_pvalue
#23 frequentist_add_info
#24 frequentist_add_beta_1
#25 frequentist_add_se_1
#26 comment

# To add STRAND and IMP
# header list obtained from the following code,
# gunzip -c /data/jampet/upload-20170920/INTERVAL_inf1_CXCL1___P09341_chr_merged.gz | \
# head -1 | \
# awk '{gsub(/ /, "\n",$0)};1'| awk '{print "#" NR, $1}'

# drop header but avoid order.awk so as to make up for addinfo.awk

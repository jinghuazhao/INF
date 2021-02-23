# export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107/
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
  }' ${HGI}/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.txt.gz_1.0E-5.txt
  #  ${HGI}/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37_1.0E-5.txt
) | \
gzip -f > ${INF}/HGI/gsmr_C2.txt.gz

(
  echo "SNP A1 A2 freq b se p N"
# gunzip -c $HGI/eur/COVID19_HGI_C2_ALL_eur_leave_23andme_20201020.b37.txt.gz | \
# gunzip -c $HGI/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
  gunzip -c $HGI/COVID19_HGI_C2_ALL_20210107.b37.txt.gz | \
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
gzip -f > ${INF}/HGI/gsmr_C2.txt.gz

#  head -1 ${HGI}/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.txt.gz_1.0E-5.txt | sed 's/\t/\n/g' | awk '{print "#" NR, $1}'
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

if [ -f ${INF}/HGI/INF1_C2.gsmr ]; then rm ${INF}/HGI/INF1_C2.gsmr; fi
(
  cat ${INF}/HGI/gsmr_C2*.gsmr | \
  head -1
  ls ${INF}/HGI/gsmr_C2*gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > ${INF}/HGI/INF1_C2.gsmr

export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20201215/results/20210107/eur
# gunzip -c $HGI/COVID19_HGI_C2_ALL_eur_20210107.txt.gz | head -1 | tr '\t' '\n' | awk '{print "#" NR,$1}'
#1 #CHR
#2 POS
#3 REF
#4 ALT
#5 SNP
#6 BQC19_EUR_AF_Allele2
#7 BQC19_EUR_AF_fc
#8 BQC19_EUR_N
#9 BelCovid_EUR_AF_Allele2
#10 BelCovid_EUR_AF_fc
#11 BelCovid_EUR_N
#12 BioVU_EUR_AF_Allele2
#13 BioVU_EUR_AF_fc
#14 BioVU_EUR_N
#15 CCPM_EUR_AF_Allele2
#16 CCPM_EUR_AF_fc
#17 CCPM_EUR_N
#18 CU_EUR_AF_Allele2
#19 CU_EUR_AF_fc
#20 CU_EUR_N
#21 DECODE_EUR_AF_Allele2
#22 DECODE_EUR_AF_fc
#23 DECODE_EUR_N
#24 EstBB_EUR_AF_Allele2
#25 EstBB_EUR_AF_fc
#26 EstBB_EUR_N
#27 FinnGen_FIN_AF_Allele2
#28 FinnGen_FIN_AF_fc
#29 FinnGen_FIN_N
#30 GCAT_EUR_AF_Allele2
#31 GCAT_EUR_AF_fc
#32 GCAT_EUR_N
#33 GENCOVID_EUR_AF_Allele2
#34 GENCOVID_EUR_AF_fc
#35 GENCOVID_EUR_N
#36 GFG_EUR_AF_Allele2
#37 GFG_EUR_AF_fc
#38 GFG_EUR_N
#39 GHS_Freeze_145_EUR_AF_Allele2
#40 GHS_Freeze_145_EUR_AF_fc
#41 GHS_Freeze_145_EUR_N
#42 Genotek_EUR_AF_Allele2
#43 Genotek_EUR_AF_fc
#44 Genotek_EUR_N
#45 INTERVAL_EUR_AF_Allele2
#46 INTERVAL_EUR_AF_fc
#47 INTERVAL_EUR_N
#48 LGDB_EUR_AF_Allele2
#49 LGDB_EUR_AF_fc
#50 LGDB_EUR_N
#51 Lifelines_EUR_AF_Allele2
#52 Lifelines_EUR_AF_fc
#53 Lifelines_EUR_N
#54 SINAI_COVID_EUR_AF_Allele2
#55 SINAI_COVID_EUR_AF_fc
#56 SINAI_COVID_EUR_N
#57 Stanford_EUR_AF_Allele2
#58 Stanford_EUR_AF_fc
#59 Stanford_EUR_N
#60 TOPMed_CHRIS10K_EUR_AF_Allele2
#61 TOPMed_CHRIS10K_EUR_AF_fc
#62 TOPMed_CHRIS10K_EUR_N
#63 TOPMed_Gardena_EUR_AF_Allele2
#64 TOPMed_Gardena_EUR_AF_fc
#65 TOPMed_Gardena_EUR_N
#66 UCLA_EUR_AF_Allele2
#67 UCLA_EUR_AF_fc
#68 UCLA_EUR_N
#69 UKBB_EUR_AF_Allele2
#70 UKBB_EUR_AF_fc
#71 UKBB_EUR_N
#72 SPGRX_EUR_AF_Allele2
#73 SPGRX_EUR_AF_fc
#74 SPGRX_EUR_N
#75 MVP_EUR_AF_Allele2
#76 MVP_EUR_AF_fc
#77 MVP_EUR_N
#78 genomicsengland100kgp_EUR_AF_Allele2
#79 genomicsengland100kgp_EUR_AF_fc
#80 genomicsengland100kgp_EUR_N
#81 Helix_EUR_AF_Allele2
#82 Helix_EUR_AF_fc
#83 Helix_EUR_N
#84 MGI_EUR_AF_Allele2
#85 MGI_EUR_AF_fc
#86 MGI_EUR_N
#87 NTR_EUR_AF_Allele2
#88 NTR_EUR_AF_fc
#89 NTR_EUR_N
#90 PHBB_EUR_AF_Allele2
#91 PHBB_EUR_AF_fc
#92 PHBB_EUR_N
#93 Ancestry_EUR_AF_Allele2
#94 Ancestry_EUR_AF_fc
#95 Ancestry_EUR_N
#96 23ANDME_EUR_AF_Allele2
#97 23ANDME_EUR_AF_fc
#98 23ANDME_EUR_N
#99 idipaz24genetics_EUR_AF_Allele2
#100 idipaz24genetics_EUR_AF_fc
#101 idipaz24genetics_EUR_N
#102 Amsterdam_UMC_COVID_study_group_EUR_AF_Allele2
#103 Amsterdam_UMC_COVID_study_group_EUR_AF_fc
#104 Amsterdam_UMC_COVID_study_group_EUR_N
#105 HOSTAGE_EUR_AF_Allele2
#106 HOSTAGE_EUR_AF_fc
#107 HOSTAGE_EUR_N
#108 SweCovid_EUR_AF_Allele2
#109 SweCovid_EUR_AF_fc
#110 SweCovid_EUR_N
#111 genomicc_EUR_AF_Allele2
#112 genomicc_EUR_AF_fc
#113 genomicc_EUR_N
#114 all_meta_N
#115 all_inv_var_meta_beta
#116 all_inv_var_meta_sebeta
#117 all_inv_var_meta_p
#118 all_inv_var_het_p
#119 all_meta_sample_N
#120 all_meta_AF
#121 rsid

(
  echo SNP A1 A2 freq b se p N
  gunzip -c $HGI/COVID19_HGI_C2_ALL_eur_20210107.txt.gz | \
  awk '
  {
    CHR=$1
    POS=$2
    a1=$4
    a2=$3
    if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
    else snpid="chr" CHR ":" POS "_" a1 "_" a2
    if (NR>1) print snpid, a1, a2, $120, $115, $116, $117, $119
  }'
) | \
gzip -f > ${INF}/HGI/gsmr_C2.txt.gz

#!/usr/bin/bash

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
  }' ${HGI}/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.txt.gz_1.0E-5.txt
) | \
gzip -f > ${INF}/HGI/gsmr_B2.txt.gz

(
  echo "SNP A1 A2 freq b se p N"
  gunzip -c $HGI/COVID19_HGI_B2_ALL_20210107.b37.txt.gz | \
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
gzip -f > ${INF}/HGI/gsmr_B2.txt.gz

if [ -f ${INF}/HGI/INF1_B2.gsmr ]; then rm ${INF}/HGI/INF1_B2.gsmr; fi
(
  cat ${INF}/HGI/gsmr_B2*.gsmr | \
  head -1
  ls ${INF}/HGI/gsmr_B2*gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > ${INF}/HGI/INF1_B2.gsmr

export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-analysis/20201215/results/20210107/eur
gunzip -c $HGI/COVID19_HGI_B2_ALL_eur_20210107.txt.gz | head -1 | tr '\t' '\n' | awk '{print "#" NR,$1}'
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
#12 CU_EUR_AF_Allele2
#13 CU_EUR_AF_fc
#14 CU_EUR_N
#15 EstBB_EUR_AF_Allele2
#16 EstBB_EUR_AF_fc
#17 EstBB_EUR_N
#18 FinnGen_FIN_AF_Allele2
#19 FinnGen_FIN_AF_fc
#20 FinnGen_FIN_N
#21 GENCOVID_EUR_AF_Allele2
#22 GENCOVID_EUR_AF_fc
#23 GENCOVID_EUR_N
#24 GHS_Freeze_145_EUR_AF_Allele2
#25 GHS_Freeze_145_EUR_AF_fc
#26 GHS_Freeze_145_EUR_N
#27 LGDB_EUR_AF_Allele2
#28 LGDB_EUR_AF_fc
#29 LGDB_EUR_N
#30 UCLA_EUR_AF_Allele2
#31 UCLA_EUR_AF_fc
#32 UCLA_EUR_N
#33 UKBB_EUR_AF_Allele2
#34 UKBB_EUR_AF_fc
#35 UKBB_EUR_N
#36 idipaz24genetics_EUR_AF_Allele2
#37 idipaz24genetics_EUR_AF_fc
#38 idipaz24genetics_EUR_N
#39 Amsterdam_UMC_COVID_study_group_EUR_AF_Allele2
#40 Amsterdam_UMC_COVID_study_group_EUR_AF_fc
#41 Amsterdam_UMC_COVID_study_group_EUR_N
#42 SPGRX_EUR_AF_Allele2
#43 SPGRX_EUR_AF_fc
#44 SPGRX_EUR_N
#45 DECODE_EUR_AF_Allele2
#46 DECODE_EUR_AF_fc
#47 DECODE_EUR_N
#48 MVP_EUR_AF_Allele2
#49 MVP_EUR_AF_fc
#50 MVP_EUR_N
#51 HOSTAGE_EUR_AF_Allele2
#52 HOSTAGE_EUR_AF_fc
#53 HOSTAGE_EUR_N
#54 23ANDME_EUR_AF_Allele2
#55 23ANDME_EUR_AF_fc
#56 23ANDME_EUR_N
#57 BoSCO_EUR_AF_Allele2
#58 BoSCO_EUR_AF_fc
#59 BoSCO_EUR_N
#60 FHoGID_EUR_AF_Allele2
#61 FHoGID_EUR_AF_fc
#62 FHoGID_EUR_N
#63 Ancestry_EUR_AF_Allele2
#64 Ancestry_EUR_AF_fc
#65 Ancestry_EUR_N
#66 SweCovid_EUR_AF_Allele2
#67 SweCovid_EUR_AF_fc
#68 SweCovid_EUR_N
#69 genomicc_EUR_AF_Allele2
#70 genomicc_EUR_AF_fc
#71 genomicc_EUR_N
#72 all_meta_N
#73 all_inv_var_meta_beta
#74 all_inv_var_meta_sebeta
#75 all_inv_var_meta_p
#76 all_inv_var_het_p
#77 all_meta_sample_N
#78 all_meta_AF
#79 rsid

(
  echo SNP A1 A2 freq b se p N
  gunzip -c $HGI/COVID19_HGI_B2_ALL_eur_20210107.txt.gz | \
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
gzip -f > ${INF}/HGI/gsmr_B2.txt.gz

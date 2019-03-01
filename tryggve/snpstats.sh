# 1-3-2019 JHZ

function EGCUT()
{
(
   awk "NR>11" /data/anekal/EGCUT_INF.snp-stats | \
   awk -vOFS="\t" '{
      if(NR==1) print "SNPID","N","MAF","HWE","info";
      else {
        CHR=$3
        POS=$4
        a1=$6
        a2=$7
        N=$26
        MAF=$14
        HWE=$8
        info=$18
        if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
        else snpid="chr" CHR ":" POS "_" a1 "_" a2
        print snpid, N, MAF, HWE, info
      }
    }'
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > snpstats/EGCUT.snpstats.gz
}

# awk 'NR>11' /data/anekal/EGCUT_INF.snp-stats  | \
# head -1 | \
# awk '{gsub(/\t/, "\n",$0);print}'| \
# awk '{print "#" NR, $1}'

#1 alternate_ids
#2 rsid
#3 chromosome
#4 position
#5 alleleA
#6 alleleB
#7 comment
#8 HW_exact_p_value
#9 HW_lrt_p_value
#10 alleleA_count
#11 alleleB_count
#12 alleleA_frequency
#13 alleleB_frequency
#14 minor_allele_frequency
#15 minor_allele
#16 major_allele
#17 info
#18 impute_info
#19 missing_proportion
#20 A
#21 B
#22 AA
#23 AB
#24 BB
#25 NULL
#26 total

function INTERVAL()
{
(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/INTERVAL/impute_{}_interval.snpstats.gz | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        POS=\$4
        a1=\$5
        a2=\$6
        N=\$9+\$10+\$11
        MAF=\$15
        HWE=\$16
        info=\$19
        if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
        else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
        print snpid, N, MAF, HWE, info
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > snpstats/INTERVAL.snpstats.gz
}

# gunzip -c impute_1_interval.snpstats.gz | \
# head -1 | \
# awk '{gsub(/\t/, "\n",$0)};1'| awk '{print "#" NR, $1}'

#1 SNPID
#2 RSID
#3 chromosome
#4 position
#5 A_allele
#6 B_allele
#7 minor_allele
#8 major_allele
#9 AA
#10 AB
#11 BB
#12 AA_calls
#13 AB_calls
#14 BB_calls
#15 MAF
#16 HWE
#17 missing
#18 missing_calls
#19 information

function ORCADES()
{
(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/ORCADES/orcades_chr_1_hrc.snp-stats.gz | \
    awk 'NR>13' | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        POS=\$4
        a1=\$6
        a2=\$7
        N=\$26
        MAF=\$14
        HWE=\$8
        info=\$18
        if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
        else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
        print snpid, N, MAF, HWE, info
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > snpstats/ORCADES.snpstats.gz
}

# gunzip -c orcades_chr_1_hrc.snp-stats.gz | \
# head -14 | \
# awk '(NR==14){gsub(/\t/, "\n",$0);print}'| \
# awk '{print "#" NR, $1}'

#1 alternate_ids
#2 rsid
#3 chromosome
#4 position
#5 alleleA
#6 alleleB
#7 comment
#8 HW_exact_p_value
#9 HW_lrt_p_value
#10 alleleA_count
#11 alleleB_count
#12 alleleA_frequency
#13 alleleB_frequency
#14 minor_allele_frequency
#15 minor_allele
#16 major_allele
#17 info
#18 impute_info
#19 missing_proportion
#20 A
#21 B
#22 AA
#23 AB
#24 BB
#25 NULL
#26 total

function VIS()
{
(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/VIS/vis_chr{}_HRC.r1-1_nomono_I4.snp-stats.gz | \
    awk 'NR>13' | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        POS=\$4
        a1=\$6
        a2=\$7
        N=\$26
        MAF=\$14
        HWE=\$8
        info=\$18
        if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
        else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
        print snpid, N, MAF, HWE, info
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > snpstats/VIS.snpstats.gz
}

# gunzip -c vis_chr1_HRC.r1-1_nomono_I4.snp-stats.gz | \
# head -14 | \
# awk '(NR==14){gsub(/\t/, "\n",$0);print}'| \
# awk '{print "#" NR, $1}'

#1 alternate_ids
#2 rsid
#3 chromosome
#4 position
#5 alleleA
#6 alleleB
#7 comment
#8 HW_exact_p_value
#9 HW_lrt_p_value
#10 alleleA_count
#11 alleleB_count
#12 alleleA_frequency
#13 alleleB_frequency
#14 minor_allele_frequency
#15 minor_allele
#16 major_allele
#17 info
#18 impute_info
#19 missing_proportion
#20 A
#21 B
#22 AA
#23 AB
#24 BB
#25 NULL
#26 total

function STABILITY()
{
(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/STABILITY/STABILITY_imputed_qc_maf_005_chr1_snp_stats.txt.gz | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        POS=\$4
        a1=\$5
        a2=\$6
        N=\$9+\$10+\$11
        MAF=\$15
        HWE=\$16
        info=\$19
        if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
        else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
        print snpid, N, MAF, HWE, info
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > snpstats/INTERVAL.snpstats.gz
}

# gunzip -c STABILITY_imputed_qc_maf_005_chr1_snp_stats.txt.gz | \
# head -1 | \
# awk '{gsub(/ /, "\n",$0);print}'| \
# awk '{print "#" NR, $1}'

#1 SNPID
#2 RSID
#3 chromosome
#4 position
#5 A_allele
#6 B_allele
#7 minor_allele
#8 major_allele
#9 AA
#10 AB
#11 BB
#12 AA_calls
#13 AB_calls
#14 BB_calls
#15 MAF
#16 HWE
#17 missing
#18 missing_calls
#19 information

for s in EGCUT INTERVAL ORCADES VIS STABILITY
do
  $s
done

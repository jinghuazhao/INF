# 14-5-2019 JHZ

# ORCADES, STABILITY, VIS and MadCam actually contains INFO
# INTERVAL and KORA has SNPTEST outputs as well as qctool -snp-stats
# ORCADES, STABILITY and VIS have chromosome as NA
# BioFinder, NSPHS, RECOMBINE and STANLEY are pending on qctool -snp-stats

(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/INTERVAL/impute_{}_interval.snpstats.gz | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"SNP\",\"a1\",\"a2\";
      else {
        CHR=\$3
        sub(/^0/,\"\",CHR)
        POS=\$4
        a1=\$5
        a2=\$6
        if (a1>a2) print \"chr\" CHR \":\" POS \"_\" a2 \"_\" a1, \$1, a2, a1;
        print \"chr\" CHR \":\" POS \"_\" a1 \"_\" a2, \$1, a1, a2
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > work/INTERVAL.rsid.gz

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

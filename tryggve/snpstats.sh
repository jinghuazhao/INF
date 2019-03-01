# 1-3-2019 JHZ

function INTERVAL()
{
(
  seq 22 | \
  parallel -j1 -C' ' '
    gunzip -c /data/jinhua/data/INTERVAL/impute_{}_interval.snpstats.gz | \
    awk "{
      OFS="\t"
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        POS=\$4
        a1=\$6
        a2=\$7
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
awk '(NR==1||$3!="chromosome")' | \
sort -k1,1 | \
gzip -f > snpstats/INTERVAL.snpstats.gz
}

INTERVAL

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

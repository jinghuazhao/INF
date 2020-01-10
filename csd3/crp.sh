# 10-1-2020 JHZ

# mg/L
gunzip -c ukb/30710_raw.gwas.imputed_v3.both_sexes.tsv.bgz | \
awk '{
   OFS="\t"
   if (NR>1)
   {
     split($1,a,":")
     CHR=a[1]
     POS=a[2]
     a1=a[3]
     a2=a[4]
     if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
     else snpid="chr" CHR ":" POS "_" a1 "_" a2
     $1=snpid
   }
   print
}' | \
sort -k1,1 | \
join - <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1) > work//30710_raw

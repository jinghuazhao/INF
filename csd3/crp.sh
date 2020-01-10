# 10-1-2020 JHZ

# mg/L
gunzip -c ukb/30710_raw.gwas.imputed_v3.both_sexes.tsv.bgz |
sed '1d' | \
awk -f tryggve/order.awk | \
sort -k1,1 | \
join - <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1) > work//30710_raw

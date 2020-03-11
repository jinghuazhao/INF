# 11-3-2020 JHZ

export bfile=${INF}/INTERVAL/cardio/INTERVAL

cut -f1-3,5,6 work/INF1.merge | \
awk 'NR>1 {gsub(/chr/,"",$1);print}' | \
parallel --env bfile=$bfile -j5 -C' ' '
plink --bfile ${bfile} --ld-snp {5} --chr {1} --from-bp {2} --to-bp {3} --r2 --ld-window-r2 0.6 --out work/{4}-{5}
'

# all variants
# awk -v chr={1} -v from={2} -v to={3} "chr==\$1 && \$4 >= from && \$4 <= to" ${bfile}.bim

# exclusion list

awk '$4>=53296855 && $4<=54500000{print $2}' ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_19.bim | \
sort > ${INF}/work/NLRP2

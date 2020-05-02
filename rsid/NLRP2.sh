awk '!($3>=5329685 && $3<=54500000){print $2}' ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_19.bim > ${INF}/work/NLRP2

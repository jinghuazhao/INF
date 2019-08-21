# 21-8-2019 JHZ

export f=INTERVAL/per_chr/interval.imputed.olink.chr_
seq 22 | parallel --env f -j1 -C' ' 'awk -f csd3/ld-map.awk ${f}{}.bim > work/INTERVAL-{}.map'
seq 22 | parallel --env f -j1 -C' ' 'qctool -g ${f}{}.bgen -map-id-data work/INTERVAL-{}.map -og work/INTERVAL-{}.bgen'
seq 22 | parallel --env f -j1 -C' ' 'bgenix -g ${f}{}.bgen -index'

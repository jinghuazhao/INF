# 22-8-2019 JHZ

export f=INTERVAL/per_chr/interval.imputed.olink.chr_
export TMPDIR=/rds/user/jhz22/hpc-work/work

seq 22 | parallel --env f -j1 -C' ' 'awk -f csd3/ld-map.awk ${f}{}.bim > work/INTERVAL-{}.map'
seq 22 | parallel --env f -j1 -C' ' 'qctool -g ${f}{}.bgen -map-id-data work/INTERVAL-{}.map -og work/INTERVAL-{}.bgen'
seq 22 | parallel --env f -j1 -C' ' 'bgenix -g work/INTERVAL-{}.bgen -index'

cat-bgen -g $(seq 22|awk '{printf "work/INTERVAL-" $1 ".bgen "}') -og work/INTERVAL.bgen

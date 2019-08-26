# 26-8-2019 JHZ

export f=INTERVAL/per_chr/interval.imputed.olink.chr_
export TMPDIR=/rds/user/jhz22/hpc-work/work

# snpid - rsid mapping

seq 22 | parallel --env f -j1 -C' ' '(cat csd3/ld-map.hdr;awk -f csd3/ld-map.awk ${f}{}.bim) > work/INTERVAL-{}.map'
cat $(seq 22 | awk '{printf "work/INTERVAL-" $1 ".map "}') | awk -vOFS="\t" '{print $8,$7}' | sort -k1,1 | gzip -f > INTERVAL/INTERVAL.rsid.gz

# bgen + bgen.bgi

seq 22 | parallel --env f -j1 -C' ' 'qctool -g ${f}{}.bgen -map-id-data work/INTERVAL-{}.map -og work/INTERVAL-{}.bgen'
seq 22 | parallel --env f -j1 -C' ' 'bgenix -g work/INTERVAL-{}.bgen -index'

cat-bgen -g $(seq 22|awk '{printf "work/INTERVAL-" $1 ".bgen "}') -og work/INTERVAL.bgen -clobber

# bed + bim + fam

# NLRP12 - chr19:54,296,855-54,311,176
qctool -g work/INTERVAL.bgen -s ${f}10.sample -excl-range 19:54296855-54500000 -ofiletype binary_ped -og work/INTERVAL

### NOTE bgen, rsid, bed+bim+fam are moved to INTERVAL directory

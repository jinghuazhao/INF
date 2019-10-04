# 4-10-2019 JHZ

export f=INTERVAL/per_chr/interval.imputed.olink.chr_
export TMPDIR=/rds/user/jhz22/hpc-work/work

# bgen in INTERVAL/ directory

## snpid - rsid mapping

seq 22 | parallel --env f -j1 -C' ' '(cat csd3/ld.hdr;awk -f csd3/ld.awk ${f}{}.bim) > INTERVAL/INTERVAL-{}.map'

## bgen

seq 22 | parallel --env f -j1 -C' ' 'qctool -g ${f}{}.bgen -map-id-data INTERVAL/INTERVAL-{}.map -og INTERVAL/INTERVAL-{}.bgen'

cat-bgen -g $(seq 22|awk '{printf "INTERVAL/INTERVAL-" $1 ".bgen "}') -og INTERVAL.bgen -clobber

# version with no duplicates and extended (chr19:53296855-54500000) NLRP12 region (chr19:54296855-54311176) in work/ directory

## duplicates

module load plink/2.00-alpha
plink2 --bgen INTERVAL/INTERVAL.bgen -sample INTERVAL/o5000-inf1-outlier_in-r2.sample --rm-dup force-first list --out dup

## snpid - rsid with no duplicates and NLRP12
## INTERVAL.rsid can alsb be from fully imputed/genotyped data (tryggve/snpid-rsid.sh, +34684915) for annotation.
## INTERVAL.snpid will be used for data extraction from INTERVAL as such.

join -v1 <(cat $(seq 22 | \
           awk '{printf "INTERVAL/INTERVAL-" $1 ".map "}') | \
           awk -vOFS="\t" '!($9==19 && $10 >= 53296855 && $10 < 54500000){print $8,$7}' | \
           sort -k1,1) \
         <(sort -k1,1 INTERVAL/dup.rmdup.list) \
         > work/INTERVAL.rsid
cut -d' ' -f1 work/INTERVAL.rsid > work/INTERVAL.snpid

## bgen -- ~35hr, infeasible under HPC
## qctool -g INTERVAL.bgen -incl-rsids work/INTERVAL.snpid -threads 5 -og work/INTERVAL.bgen

## bed + bim + fam

## The bgen file below has no duplicates
### plink2 --bgen work/INTERVAL.bgen -sample INTERVAL/o5000-inf1-outlier_in-r2.sample --make-bed --out INTERVAL/INTERVAL

qctool -g work/INTERVAL.bgen -s INTERVAL/o5000-inf1-outlier_in-r2.sample -threads 5 -ofiletype binary_ped -og INTERVAL/INTERVAL

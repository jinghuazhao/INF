# 22-4-2020 JHZ22

function sentinels()
# f10 is from S10
# f11 is S11
{
  join -12 -t$'\t' <(sed '1d' f10 | cut -f1,4 | sort -k2,2) <(sed '1d' f11 | cut -f1-3 | sort -k1,1) | \
  sort -k2,2 | awk -v FS="\t" -vOFS="\t" '{split($3,a,":");$3 ="chr" a[1] ":" a[2] "_" a[3]};1' | \
  join -12 -t$'\t' - <(sed '1d' doc/olink.inf.panel.annot.tsv | cut -f1-3 | \
                       awk -vFS="\t" -vOFS="\t" '{gsub(/\"/,"",$3);gsub(/\"/,"",$2);print $3,$2,$1}' | \
                       sort -k1,1) | \
  join -t$'\t' -12 <(sort -k2,2 work/inf1.tmp) - | \
  sort -k2,2 | \
  join -12 -t$'\t' - <(sed '1d' work/INF1.merge | cut -f5,6 | sort -k1,1)
}

sentinels | \
awk -vFS="\t" -vOFS="\t" '$4==$8{print $5,$2,$1,$7," ","Folkersen, et al (2020)"}' | \
xsel -i

module load plink/2.00-alpha
(
  sentinels | \
  awk -vFS="\t" -vOFS="\t" '$4!=$8{print $8,$1,$2,$3,$4,$5,$6,$7}' | \
  sort -k1,1 | \
  join <(sort -k1,1 work/INF1.merge.rsid | awk -vOFS="\t" '{split($1,a,":");sub(/chr/,"",a[1]);print $1,$2,a[1]}') - | \
  parallel -C' ' '
    echo {1} {2} {3} {4} {5} {6} {7} {8} {9}
    plink2 --bfile INTERVAL/per_chr/interval.imputed.olink.chr_{3} --ld {2} {8}
  '
) > cvd1.log

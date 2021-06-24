#!/usr/bin/bash

export chr=1
export start=159175353
export end=159525679
export M=1e6

# SCALLOP/INF
for prot in MCP.1 MCP.2 MCP.3 MCP.4
do
  gunzip -c METAL/${prot}-1.tbl.gz | \
  awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
  {
    if ($1 == chr && $2 >= start - M && $2 <= end + M) 
    {
      split($3,a,"_")
      print a[1],$1,$2,$10,$3,toupper($4),toupper($5)
    }
  }' | \
  sort -k1,1 | \
  join -12 -21 work/snp_pos - | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > work/${prot}-pQTL.beta
done

# Monoocyte count
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/mono.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > work/mono-QTL.beta

join -j5 <(sort -k5,5 work/MCP.1-pQTL.beta) <(sort -k5,5 work/MCP.2-pQTL.beta) | \
join -25 - <(sort -k5,5 work/MCP.3-pQTL.beta) | \
join -25 - <(sort -k5,5 work/MCP.4-pQTL.beta) | \
join -25 - <(sort -k5,5 work/mono-QTL.beta) | \
awk -vOFS="\t" '
{
  if($6!=12) $11=-$11
  if($6!=18) $17=-$17
  if($6!=24) $23=-$23
  if($6!=30) $29=-$29
  print $1,$2,$3,$4,$5,$11,$17,$23,$29
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n > work/rs12075-beta.gassoc

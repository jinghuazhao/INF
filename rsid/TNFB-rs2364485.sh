#!/usr/bin/bash

export chr=12
export start=6499532
export end=6519837
export M=1e6
export gene=LTBR
export prot=TNFB

echo Multiple sclerosis
# rs2364485 and rs1800693
zgrep -v -w -e NaN -e Infinity -e 0.0 annotate/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz | \
sed '1d' | \
awk '!/^$/' | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($1==chr && $2>=start-M && $2<=end+M)
  {
    if ($4<$5) snpid="chr" $1 ":" $2 "_" $4 "_" $5;
    else snpid="chr" $1 ":" $2 "_" $5 "_" $4
    print $3,$1,$2,$7/$8,snpid,$5,$4
  }
}' > work/${prot}-QTL.lz

echo SCALLOP/INF
gunzip -c METAL/${prot}-1.tbl.gz | \
sed '1d' | \
awk -vFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
{
  if ($1 == chr && $2 >= start-M && $2 <= end+M)
  {
    split($3,a,"_")
    print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5)
  }
}' | \

sort -k1,1 | \
join -12 -21 work/snp_pos - | \
awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > work/${prot}-pQTL.lz

echo eQTLGen
zgrep -w ${gene} eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > work/${prot}-eQTL.lz

join -j5 <(sort -k5,5 work/${prot}-QTL.lz) <(sort -k5,5 work/${prot}-pQTL.lz) | \
join -25 - <(sort -k5,5 work/${prot}-eQTL.lz) | \
awk -vOFS="\t" '
{
  if($6!=12) $11=-$11
  if($6!=18) $17=-$17
  print $1,$2,$3,$4,$5,$11,$17
}' | \
awk 'a[$1]++==0' | \
sort -k3,3n -k4,4n > work/${prot}.gassoc

cut -f1 work/${prot}.gassoc > work/${prot}.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/${prot}.snpid --r square --out work/${prot}

R --no-save -q <<END
  prot <- Sys.getenv("prot")
  d <- read.table(paste0(file.path("work",prot),".gassoc"),col.names=c("snpid","marker","chr","pos","QTL","pQTL","eQTL"))
  markers <- d[c("marker","chr","pos")]
  z <- d[c("QTL","pQTL","eQTL")]
  rownames(z) <- with(d,marker)
  ld <- read.table(paste0(file.path("work",prot),".ld"),col.names=with(d,marker),row.names=with(d,marker))
  library(gassocplot)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("QTL","pQTL","eQTL"), ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-rs385076.png"), 3, width=8, height=13)
END

# Multiple_sclerosis
chrom	pos	rsid	other_allele	effect_allele	p	beta	se	OR	OR_lower	OR_upper
1	1118275	rs61733845	G	A	0.762	-0.011971371781219958	0.03952833128852883	0.9881	0.91443681941827404	1.0676971763025762

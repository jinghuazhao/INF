#!/usr/bin/bash

export chr=12
export start=6400000
export end=6519837
export M=0
export gene=LTBR
export prot=TNFB

echo Multiple sclerosis
# rs1800693 chr12:6440009
# rs2364485 chr12:6514963
# 2018, the rsid is incompleete
zgrep -w -e 6440009 -e 6514963 data/discovery_metav3.0.meta.gz
# 12 6440009 rs1800693 T C 14 1.017e-13 0.8808
# 12 6514963 rs2364485 C A 15 5.778e-06 0.9041
gunzip -c data/discovery_metav3.0.meta.gz | \
sed '1d' | \
awk '!/^$/' | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($1==chr && $2>=start-M && $2<=end+M)
  {
    if ($4<$5) snpid="chr" $1 ":" $2 "_" $4 "_" $5;
    else snpid="chr" $1 ":" $2 "_" $5 "_" $4
    print "chr"$1":"$2,$1,$2,$7,snpid,$4,$5
  }
}' | \
sort -k1,1 | \
join -12 -21 work/snp_pos - | \
awk -vOFS="\t" '{print $2,$3,$4,$5,$6,$7,$8}' > work/${prot}-QTL.lz

echo 2013
zgrep -w -e rs1800693 -e rs2364485 data/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz
zgrep -v -w -e NaN -e Infinity -e 0.0 data/imsgc_2013_24076602_ms_efo0003885_1_ichip.sumstats.tsv.gz | \
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
}' > work/${prot}-2013.lz

echo SCALLOP/INF
# chr12:6514963_A_C
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
zgrep -w ${gene} data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > work/${prot}-eQTL.lz

join -a2 -e "NA" -o2.5,2.1,2.2,2.3,1.4,1.6,1.7,2.1,2.2,2.3,2.4,2.6,2.7 \
     -j5 <(sort -k5,5 work/${prot}-QTL.lz) <(sort -k5,5 work/${prot}-pQTL.lz) | \
join -a1 -e "NA" -25 -o1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.1,2.2,2.3,2.4,2.6,2.7 \
     - <(sort -k5,5 work/${prot}-eQTL.lz) | \
awk -vOFS="\t" '
{
# 2013 with beta/se
# if($6!="NA" && $12!="NA" && $11!="NA" && $6!=12) $11=-$11
# if($6!="NA" && $18!="NA" && $17!="NA" && $6!=18) $17=-$17
  print $1,$2,$3,$4,$5,$11,$17
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n > work/${prot}.z

cut -f1 work/${prot}.gassoc > work/${prot}.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/${prot}.snpid --make-bed --out work/${prot}
cut -f2 work/${prot}.bim > work/${prot}.snpid
plink --bfile work/${prot} --extract work/${prot}.snpid --r square --out work/${prot}
grep -f work/${prot}.snpid work/${prot}.z > work/${prot}.gassoc

R --no-save -q <<END
  prot <- Sys.getenv("prot")
  d <- read.table(paste0(file.path("work",prot),".gassoc"),
                         col.names=c("snpid","marker","chr","pos","Multiple_sclerosis","pQTL","eQTL"))
  d <- within(d,{Multiple_sclerosis <- qnorm(Multiple_sclerosis/2)})
  markers <- d[c("marker","chr","pos")]
  z <- d[c("Multiple_sclerosis","pQTL","eQTL")]
  rownames(z) <- with(d,marker)
  ld <- read.table(paste0(file.path("work",prot),".ld"),col.names=with(d,marker),row.names=with(d,marker))
  library(gassocplot)
  for (rsid in c("rs1800693","rs2364485"))
  {
    sap <- stack_assoc_plot(markers, z, ld, top.marker=rsid, traits = c("Multiple_sclerosis","pQTL","eQTL"), ylab = "-log10(P)", legend=TRUE)
    stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-",rsid,".png"), 3, width=8, height=13)
  }
  z <- within(z,{Multiple_sclerosis <- NA})
# rs1800693: P 2e-47 Z -14.46555
  pos1 <- with(d,marker=="rs1800693")
  z[pos1,"Multiple_sclerosis"] <- -14.46555
# rs2364485: P 2e-20 Z  -9.26234
  pos2 <- with(d,marker=="rs2364485")
  z[pos2,"Multiple_sclerosis"] <- -9.26234
  rsid <- "rs2364485"
  sap <- stack_assoc_plot(markers, z, ld, top.marker=rsid, traits = c("Multiple_sclerosis","pQTL","eQTL"), ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-",rsid,"-fixed.png"), 3, width=8, height=13)
END

# Multiple_sclerosis

# 2018
CHR BP SNP A1 A2 N P OR
1 11154 chr1:11154 C A 4 0.7911 0.9818

# 2013
# chrom	pos	rsid	other_allele	effect_allele	p	beta	se	OR	OR_lower	OR_upper
# 12      6440009 rs1800693       A       G       6.92e-16        0.13453089295760606     0.01666652509712173     1.144   1.1072334356041962      1.1819874273267836

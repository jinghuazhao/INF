#!/usr/bin/bash

export chr=9
export start=136155000
export end=136155000
export M=1e6

# SCALLOP/INF
for prot in CCL25 CX3CL1 LIF.R
do
  gunzip -c METAL/${prot}-1.tbl.gz | \
  awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
  {
    if ($1 == chr && $2 >= start - M && $2 <= end + M) 
    {
      split($3,a,"_")
      print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5),$11
    }
  }' | \
  sort -k1,1 | \
  join -12 -21 work/snp_pos - | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8, $9}' > work/HGI/${prot}-pQTL.lz
done

# Covid-19
export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20200915/results/20201020
module load gcc/6
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($3<$4) snpid="chr" $1 ":" $2 "_" $3 "_" $4;
  else snpid="chr" $1 ":" $2 "_" $4 "_" $3
  if($1==chr && $2>=start-M && $2 <=end+M) print $13,$1,$2,$7/$8,snpid,$3,$4,$8
}' ${HGI}/COVID19_HGI_C2_ALL_20201020.b37_1.0E-5.txt > work/HGI/HGI-QTL.lz

join -j5 <(sort -k5,5 work/HGI/CCL25-pQTL.lz) <(sort -k5,5 work/HGI/CX3CL1-pQTL.lz) | \
join -25 - <(sort -k5,5 work/HGI/LIF.R-pQTL.lz) | \
join -25 - <(sort -k5,5 work/HGI/HGI-QTL.lz) | \
awk -vOFS="\t" '
{
  if($6!=$13) $12=-$12
  if($6!=$20) $19=-$19
  if($6!=$27) $26=-$26
  print $1,$2,$3,$4,$5,$12,$19,$26,$8,$15,$22,$29
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n | \
grep -v chr9:136156051_A_G > work/HGI/rs635634.gassoc

cut -f1 work/HGI/rs635634.gassoc | \
grep -v chr9:136156051_A_G > work/HGI/rs635634.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/HGI/rs635634.snpid --r square --out work/HGI/rs635634

R --no-save -q <<END
  library(gassocplot)
  c3 <- c("CCL25","CX3CL1","LIF.R")
  c4 <- c(c3,"C2")
  d <- read.table("work/HGI/rs635634.gassoc",col.names=c("snpid","marker","chr","pos",c3,"C2",paste0("se_",c4)))
  markernames <- with(d,as.character(marker))
  markers <- d[c("marker","chr","pos")]
  z <- d[c4]
  rownames(z) <- markernames
  ld <- read.table("work/HGI/rs635634.ld",col.names=markernames,row.names=markernames)
  sap <- stack_assoc_plot(markers, z, ld, traits = c4, ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, "work/HGI/rs635634.png", 5, width=8, height=13, dpi=300)
  library(MendelianRandomization)
  by <- d[["C2"]]
  byse <- d[["se_C2"]]
  for (p in c("CCL25","CX3CL1", "LIF.R"))
  {
    bx <- d[[p]]
    bxse <- d[[paste0("se_",p)]]
    bse <- mr_input(bx,bxse,by,byse,snps=markernames)
    mr_plot(bse)
    print(mr_ivw(bse))
    print(mr_egger(bse))
  }
  bx <- cbind(CCL25=d[["CCL25"]],CX3CL1=d[["CX3CL1"]],LIF.R=d[["LIF.R"]])
  bxse <- cbind(se_CCL25=d[["se_CCL25"]],se_CX3CL1=d[["se_CX3CL1"]],se_LIF.R=d[["se_LIF.R"]])
  mvbse <- mr_mvinput(bx,bxse,by,byse,snps=markernames)
  mr_mvivw(mvbse)
  mr_mvegger(mvbse)
  mr_mvmedian(mvbse)
  mr_plot(mvbse)
END

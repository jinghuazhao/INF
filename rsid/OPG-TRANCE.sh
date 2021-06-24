#!/usr/bin/bash

export chr=8
export start=119900000
export end=120300000
export M=0

echo SCALLOP/INF
# chr8:120081031_C_T rs2247769
# chr8:120201029_C_T rs2468187
for prot in OPG TRANCE
do
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
  awk 'a[$6]++==0' | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > work/${prot}-pQTL.lz
done

export prot=OPG-TRANCE
join -a2 -e "NA" -o2.5,2.1,2.2,2.3,1.4,1.6,1.7,2.1,2.2,2.3,2.4,2.6,2.7 \
     -j5 <(sort -k5,5 work/OPG-pQTL.lz) <(sort -k5,5 work/TRANCE-pQTL.lz) | \
awk -vOFS="\t" '
{
  if($6!="NA" && $12!="NA" && $11!="NA" && $6!=12) $11=-$11
  print $1,$2,$3,$4,$5,$11
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n > work/${prot}.z

cut -f1 work/${prot}.z > work/${prot}.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/${prot}.snpid --make-bed --out work/${prot}
cut -f2 work/${prot}.bim > work/${prot}.snpid
plink --bfile work/${prot} --extract work/${prot}.snpid --r square --out work/${prot}
grep -f work/${prot}.snpid work/${prot}.z > work/${prot}.gassoc

R --no-save -q <<END
  prot <- Sys.getenv("prot")
  d <- read.table(paste0(file.path("work",prot),".gassoc"),
                         col.names=c("snpid","marker","chr","pos","OPG","TRANCE"))
  markers <- d[c("marker","chr","pos")]
  z <- d[c("OPG","TRANCE")]
  rownames(z) <- with(d,marker)
  ld <- read.table(paste0(file.path("work",prot),".ld"),col.names=with(d,marker),row.names=with(d,marker))
  library(gassocplot)
  for (rsid in c("rs2247769","rs2468187"))
  {
    sap <- stack_assoc_plot(markers, z, ld, top.marker=rsid, traits = c("OPG","TRANCE"), ylab = "-log10(P)", legend=TRUE)
    stack_assoc_plot_save(sap, paste0(file.path("work",prot),"-",rsid,".png"), 3, width=8, height=13, dpi=300)
  }
END

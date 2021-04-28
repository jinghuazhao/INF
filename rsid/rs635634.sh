#!/usr/bin/bash

export chr=9
export start=136155000
export end=136155000
export M=1e6

function lz()
{
# SCALLOP/INF
for prot in CCL25 CX3CL1 LIF.R
do
  gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
  awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
  {
    if ($1 == chr && $2 >= start - M && $2 <= end + M) 
    {
      split($3,a,"_")
      print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5),$11
    }
  }' | \
  sort -k1,1 | \
  join -12 -21 ${INF}/work/snp_pos - | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8, $9}' > ${INF}/HGI/${prot}-pQTL.lz
done

# Covid-19
export HGI=~/rds/results/public/gwas/covid19/hgi/covid19-hg-public/20201215/results/20210107
for trait in A2 B2 C2
do
   gunzip -c ${HGI}/COVID19_HGI_${trait}_ALL_eur_leave_23andme_20210107.b37.txt.gz | \
   awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
   {
     if ($3<$4) snpid="chr" $1 ":" $2 "_" $3 "_" $4;
     else snpid="chr" $1 ":" $2 "_" $4 "_" $3
     if($1==chr && $2>=start-M && $2 <=end+M) print $13,$1,$2,$7/$8,snpid,$3,$4,$8
   }' > ${INF}/HGI/HGI-${trait}.lz
done
}

function ld()
{
for trait in A2 B2 C2
do
  join -j5 <(sort -k5,5 ${INF}/HGI/CCL25-pQTL.lz) <(sort -k5,5 ${INF}/HGI/CX3CL1-pQTL.lz) | \
  join -25 - <(sort -k5,5 ${INF}/HGI/LIF.R-pQTL.lz) | \
  join -25 - <(sort -k5,5 ${INF}/HGI/HGI-${trait}.lz) | \
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
  grep -v chr9:136156051_A_G > ${INF}/HGI/rs635634-${trait}.gassoc

  cut -f1 ${INF}/HGI/rs635634-${trait}.gassoc | \
  grep -v chr9:136156051_A_G > ${INF}/HGI/rs635634-${trait}.snpid
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/HGI/rs635634-${trait}.snpid --make-bed --out ${INF}/HGI/rs635634-${trait}
  cut -f2 ${INF}/HGI/rs635634-${trait}.bim > ${INF}/HGI/rs635634-${trait}.left
  plink --bfile ${INF}/HGI/rs635634-${trait} --r square --out ${INF}/HGI/rs635634-${trait}
done
}

module load gcc/6

function mvmr()
{
for trait in A2 B2 C2; do export trait=${trait}
R --no-save -q <<END
  trait <- Sys.getenv("trait")
  cat(trait,"\n")
  INF <- Sys.getenv("INF")
  library(gassocplot)
  c3 <- c("CCL25","CX3CL1","LIF.R")
  c4 <- c(c3,trait)
  d <- read.table(file.path(INF,"HGI",paste0("rs635634-",trait,".gassoc")),
                  col.names=c("snpid","marker","chr","pos",c3,trait,paste0("se_",c4)))
  d <- subset(d, snpid%in% read.table(file.path(INF,"HGI",paste0("rs635634-",trait,".left")),col.names="snpid")[["snpid"]])
  markernames <- with(d,as.character(marker))
  markers <- d[c("marker","chr","pos")]
  z <- d[c4]
  rownames(z) <- markernames
  ld <- read.table(file.path(INF,"HGI",paste0("rs635634-",trait,".ld")),col.names=markernames,row.names=markernames)
  sap <- stack_assoc_plot(markers, z, ld, traits = c4, ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, file.path(INF,"HGI",paste0("rs635634-",trait,".png")), 5, width=8, height=13, dpi=300)
  library(MendelianRandomization)
  by <- d[[trait]]
  byse <- d[[paste0("se_",trait)]]
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
done
}

function ma()
(
  gunzip -c ${INF}/HGI/gsmr_C2.txt.gz | head -1
  gunzip -c ${INF}/HGI/gsmr_C2.txt.gz | sed '1d' | sort -k1,1 | join <(gunzip -c ~/SUMSTATS/snp150.snpid_rsid.gz) - | \
  cut -d' ' -f1 --complement
) > ${INF}/HGI/gcta_C2.ma

mvmr > ${INF}/HGI/rs635634-A2-B2-C2.out

R --no-save <<END
  library(ieugwasr)
  # GET
  api_query("associations/prot-a-1737/rs635634")
  # POST
  api_query("associations", query=list(rsid="rs635634", id="prot-a-1737"))
  associations(variants="rs635634", id=c("prot-a-1737"))
  rs635634 <- subset(data.frame(phewas(variants="rs635634", pval=5e-8)),grepl("prot",id))
  library(pQTLtools)
  d <- subset(merge(SomaLogic160410,rs635634,by.x="TargetFullName",by.y="trait"),entGene=="LIFR")
  d
  gwasinfo(with(d,id))
  data.frame(ieugwasr::afl2_rsid(c("rs635634")))
END

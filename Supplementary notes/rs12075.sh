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
      print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5)
    }
  }' | \
  sort -k1,1 | \
  join -12 -21 work/snp_pos - | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > work/${prot}-pQTL.lz
done

# Baseophil, monocyte, WBC counts
for trait in baso mono wbc
do
  gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/${trait}.tsv.gz | \
  awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
  {
    if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
    else snpid="chr" $3 ":" $4 "_" $6 "_" $5
    if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7/$8,snpid,$5,$6
  }' > work/${trait}-QTL.lz
done

join -j5 <(sort -k5,5 work/MCP.1-pQTL.lz) <(sort -k5,5 work/MCP.2-pQTL.lz) | \
join -25 - <(sort -k5,5 work/MCP.3-pQTL.lz) | \
join -25 - <(sort -k5,5 work/MCP.4-pQTL.lz) | \
join -25 - <(sort -k5,5 work/mono-QTL.lz) | \
awk -vOFS="\t" '
{
  if($6!=$12) $11=-$11
  if($6!=$18) $17=-$17
  if($6!=$24) $23=-$23
  if($6!=$30) $29=-$29
  print $1,$2,$3,$4,$5,$11,$17,$23,$29
}' | \
awk 'a[$1]++==0' | \
awk 'a[$2]++==0' | \
sort -k3,3n -k4,4n > work/rs12075.gassoc

cut -f1 work/rs12075.gassoc > work/rs12075.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/rs12075.snpid --r square --out work/rs12075

R --no-save -q <<END
  library(gassocplot)
  d <- read.table("work/rs12075.gassoc",col.names=c("snpid","marker","chr","pos","MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count"))
  markers <- d[c("marker","chr","pos")]
  ld <- read.table("work/rs12075.ld",col.names=with(d,marker),row.names=with(d,marker))
  z <- d[c("MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count")]
  names(z) <- names(d)[5:9] <- c(paste0("MCP-",1:4),"Monocyte count")
  rownames(z) <- with(d,marker)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("MCP-1","MCP-2","MCP-3","MCP-4","Monocyte count"), ylab = "-log10(P)",
                          x.min=159000000,x.max=159600000,legend=TRUE)
  INF <- Sys.getenv("INF")
  pdf(file.path(INF,"plots","rs12075.pdf"),width=8,height=13)
  grid::grid.draw(sap)
  dev.off()
# stack_assoc_plot_save(sap, "rs12075.png", 5, width=8, height=13, dpi=300)
END

function rs12075_beta()
## rs12075-beta.sh
{
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
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > hyprcoloc/${prot}-pQTL.beta
done

# Monoocyte count
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/mono.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > hyprcoloc/mono-QTL.beta

join -j5 <(sort -k5,5 hpprcoloc/MCP.1-pQTL.beta) <(sort -k5,5 hyprcoloc/MCP.2-pQTL.beta) | \
join -25 - <(sort -k5,5 hyprcoloc/MCP.3-pQTL.beta) | \
join -25 - <(sort -k5,5 hyprcoloc/MCP.4-pQTL.beta) | \
join -25 - <(sort -k5,5 hyprcoloc/mono-QTL.beta) | \
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
sort -k3,3n -k4,4n > hyprcoloc/rs12075-beta.gassoc
}

function rs12075_beta()
## rs12075-se.sh
{
# SCALLOP/INF
for prot in MCP.1 MCP.2 MCP.3 MCP.4
do
  gunzip -c METAL/${prot}-1.tbl.gz | \
  awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
  {
    if ($1 == chr && $2 >= start - M && $2 <= end + M)
    {
      split($3,a,"_")
      print a[1],$1,$2,$11,$3,toupper($4),toupper($5)
    }
  }' | \
  sort -k1,1 | \
  join -12 -21 work/snp_pos - | \
  awk -vOFS="\t" '{print $2, $3, $4, $5, $6, $7, $8}' > hyprcoloc/${prot}-pQTL.se
done

# Monoocyte count
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/mono.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$8,snpid,$5,$6
}' > hyprcoloc/mono-QTL.se

join -j5 <(sort -k5,5 work/MCP.1-pQTL.se) <(sort -k5,5 hyprcoloc/MCP.2-pQTL.se) | \
join -25 - <(sort -k5,5 hyprcoloc/MCP.3-pQTL.se) | \
join -25 - <(sort -k5,5 hyprcoloc/MCP.4-pQTL.se) | \
join -25 - <(sort -k5,5 hyprcoloc/mono-QTL.se) | \
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
sort -k3,3n -k4,4n > hyprcoloc/rs12075-se.gassoc
}

function hyprcoloc()
{
export r=$(awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} 'BEGIN{print chr":"start-M"-"end-M}')

if [ ! -d ${INF}/hyprcoloc ]; then mkdir ${INF}/hyprcoloc; fi
# Monoocyte count
# VARIANT ID_dbSNP49      CHR     BP      REF     ALT     EFFECT_INT      SE_INT  MLOG10P_INT     ALT_FREQ_INT    INFO_INT
export trait=mono
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/${trait}.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if(NR==1) print "chr","pos","snpid","a1","a2","af","b","se","p","n";
  else if($3==chr) print $3,$4,$2,$6,$5,$10,$7,$8,$9,173480
}' | \
bgzip -f > ${trait}.tsv.gz
module load jdk-8u141-b15-gcc-5.4.0-p4aaopt
module load gatk
module load python/3.7
cd ${HPC_WORK}/gwas2vcf
source env/bin/activate
if [ ! -f ${INF}/${trait}.vcf.gz ]; then
   python main.py --out ${INF}/${trait}.vcf.gz --data ${INF}/${trait}.tsv.gz \
                  --id ${trait} --ref human_g1k_v37.fasta --json ${INF}/rsid/gwas2vcf.json
   tabix -f ${INF}/hyprcoloc/${trait}.vcf.gz
fi
deactivate
cd -

cd ${INF}/hyprcoloc
bcftools merge -r ${r} -O z -o MCP-rs12075.vcf.gz \
         METAL/gwas2vcf/MCP.1.vcf.gz METAL/gwas2vcf/MCP.2.vcf.gz METAL/gwas2vcf/MCP.3.vcf.gz METAL/gwas2vcf/MCP.4.vcf.gz ${trait}.vcf.gz
bcftools query -f "%CHROM:%POS\_%REF\_%ALT[\t%ES][\t%SE]\n" hyprcoloc/MCP-rs12075.vcf.gz > MCP-rs12075.txt

  Rscript -e '
  n.traits <- 5
  rs12075 <- within(read.table("MCP-rs12075.txt",col.names=c("rsid",paste0("v",1:(2*n.traits)))),
             {
                 v1 <- as.numeric(v1); v2 <- as.numeric(v2)
                 v3 <- as.numeric(v3); v4 <- as.numeric(v4)
                 v5 <- as.numeric(v5); v6 <- as.numeric(v6)
                 v7 <- as.numeric(v7); v8 <- as.numeric(v8)
                 v9 <- as.numeric(v9); v10 <- as.numeric(v10)
             })
  na <- apply(is.na(rs12075),1,any)
  rs12075 <- subset(rs12075,!na)
  rsid <- with(rs12075,rsid)
  betas <- as.matrix(rs12075[paste0("v",1:n.traits)])
  ses <- as.matrix(rs12075[paste0("v",(n.traits+1):(2*n.traits))])
  colnames(betas) <- colnames(ses) <- traits <- paste0("T",1:n.traits)
  rownames(betas) <- rownames(ses) <- rsid
  hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
  traits <- c("MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count")
  d <- read.table("rs12075-beta.gassoc",col.names=c("snpid","marker","chr","pos",traits))
  markers <- d[c("marker","chr","pos")]
  betas <- as.matrix(d[traits])
  rownames(betas) <- with(d,marker)
  d <- read.table("rs12075-se.gassoc",col.names=c("snpid","marker","chr","pos",traits))
  ses <- as.matrix(d[traits])
  rownames(ses) <- with(d,marker)
  hyprcoloc::hyprcoloc(betas, ses, trait.names=colnames(ses), snp.id=with(markers,marker))
  '
  cd -
}

#Results:
#  iteration                                     traits posterior_prob
#1         1 MCP.1, MCP.2, MCP.3, MCP.4, Monocyte_count          0.978
#  regional_prob candidate_snp posterior_explained_by_snp dropped_trait
#1        0.9964       rs12075                          1            NA

zcat Whole_Blood.allpairs.txt.gz | cut -f1 | awk 'NR>1{print substr($1,1,15)}' | sort | uniq > work/wb
R --no-save -q <<END
   wb <- scan("work/wb", what="")
   g <- within(subset(grex::grex(wb),!is.na(uniprot_id)),{uniprot=trimws(uniprot_id)})
   m <- merge(gap::inf1,g,by="uniprot",all.x=TRUE)
END
}

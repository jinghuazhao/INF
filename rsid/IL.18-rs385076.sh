#!/usr/bin/bash

export chr=2
export start=3e7
export end=3.45e7
export M=1e6
export gene=IL18
export prot=IL.18

# SCALLOP/INF
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

# eQTLGen
zgrep NLRC4 eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7,snpid,$5,$6
}' > work/${prot}-eQTL.lz

# Lymphocyte count
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/lymph_p.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if ($5<$6) snpid="chr" $3 ":" $4 "_" $5 "_" $6;
  else snpid="chr" $3 ":" $4 "_" $6 "_" $5
  if($3==chr && $4>=start-M && $4 <=end+M) print $2,$3,$4,$7/$8,snpid,$5,$6
}' > work/${prot}-QTL.lz

join -j5 <(sort -k5,5 work/${prot}-pQTL.lz) <(sort -k5,5 work/${prot}-eQTL.lz) | \
join -25 - <(sort -k5,5 work/${prot}-QTL.lz) | \
awk -vOFS="\t" '
{
  if($6!=12) $11=-$11
  if($6!=18) $17=-$17
  print $1,$2,$3,$4,$5,$11,$17
}' | \
awk 'a[$1]++==0' | \
sort -k3,3n -k4,4n > work/IL.18.gassoc

cut -f1 work/IL.18.gassoc > work/IL.18.snpid
plink --bfile INTERVAL/cardio/INTERVAL --extract work/IL.18.snpid --r square --out work/IL.18

R --no-save -q <<END
  library(gassocplot)
  d <- read.table("work/IL.18.gassoc",col.names=c("snpid","marker","chr","pos","pQTL","eQTL","QTL"))
  markers <- d[c("marker","chr","pos")]
  z <- d[c("pQTL","eQTL","QTL")]
  rownames(z) <- with(d,marker)
  corr <- read.table("work/IL.18.ld",col.names=with(d,marker),row.names=with(d,marker))
  sap <- stack_assoc_plot(markers, z, corr, traits = c("pQTL","eQTL","QTL"), ylab = "-log10(P)", legend=TRUE)
  stack_assoc_plot_save(sap, "work/IL.18.png", 3, width=8, height=15)
END

# pQTL
# Chromosome	Position	MarkerName	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Effect	StdErr	log(P)	Direction	HetISq	HetChiSq	HetDf	logHetP	N
# 4	165506706	chr4:165506706_C_T	t	c	0.0099	0.0029	0.0029	0.0123	0.0257	0.0763	-0.13	-+?+pn?-??-	65.6	17.419	6	-2.105	13165

# eQTLGen
# Pvalue	SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR	BonferroniP
# 3.2717E-310	rs12230244	12	10117369	T	A	200.7534	ENSG00000172322	CLEC12A	12	10126104	34	30596	0.0	4.1662E-302

# Lymphocyte count
VARIANT	ID_dbSNP49	CHR	BP	REF	ALT	EFFECT_INT	SE_INT	MLOG10P_INT	ALT_FREQ_INT	INFO_INT
1:10177_A_AC	rs367896724	1	10177	A	AC	4.71909e-03	8.50031e-03	2.374861e-01	3.759800e-01	7.35090e-01

#!/usr/bin/bash

function INTERVAL()
{
# init -- actually two versions of RNASeq results below gives the same LTBR.lz
  export rnaseq=tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv
  export rnaseq=tensorqtl_allSNPs_MAF0.005_merged_annotated.csv
  grep -w ${rsid} ${rnaseq}
  zgrep ENSG00000256433 ${INF}/work/ensGtp.txt.gz | \
  cut -f2 | \
  zgrep -f - ${INF}/work/ensemblToGeneName.txt.gz
# LocusZoom plot
  read chr start end < st.tmp
  awk -vFS="," -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1 || ($5==chr && $6>=start && $6<=end && index($0,gene)>0)' ${rnaseq} | \
  tr "," "\t" > LTBR.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal LTBR.lz --delim tab title="INTERVAL-LTBR" \
            --markercol variant_id --pvalcol pval --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=INTERVAL --rundir .
  mv INTERVAL_chr${chr}_${bracket}.pdf INTERVAL-LTBR-cis.pdf
}

function eQTLGen()
# https://www.eqtlgen.org/trans-eqtls.html
# https://www.eqtlgen.org/cis-eqtls.html
{
  export cis=2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  export trans=2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
  read chr start end < st.tmp
  (
    gunzip -c ${cis} | \
    head -1
    gunzip -c ${cis} | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} '$3==chr && $4>=start && $4<=end && index($0,gene)>0' | \
    sort -k3,3n -k4,4n
  ) > eQTLGen.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal eQTLGen.lz --delim tab title="eQTLGen-LTBR" \
            --markercol SNP --pvalcol Pvalue --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=eQTLGen --rundir .
  mv eQTLGen_chr${chr}_${bracket}.pdf eQTLGen-LTBR-cis.pdf
}

function SCALLOP()
{
  read chr start end < st.tmp
  (
    gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
    head -1 | \
    awk '{$1=$3 " " $1;$3="";$12="P"};1' | \
    awk '{$1=$1};1'
    gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
    awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} '$1 == chr && $2 >= start && $2 <= end {$12=10^$12;print}' | \
    sort -k3,3 | \
    join -23 <(awk -vchrom=chr${chr} 'index($0,chrom)>0' INTERVAL.rsid) - | \
    sort -k2,2n -k3,3n | \
    cut -d' ' -f1 --complement
  ) | \
  tr ' ' '\t' > TNFB.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal TNFB.lz --delim tab title="SCALLOP-TNFB" \
            --markercol MarkerName --pvalcol P --chr ${chr} --start ${b1} --end ${b2} --gwas-cat whole-cat_significant-only \
            --no-date --plotonly --prefix=TNFB --rundir .
  mv TNFB_chr${chr}_${bracket}.pdf SCALLOP-TNFB-cis.pdf
}

cd work
export chr=12
export pos=6514963
export gene=LTBR
export rsid=rs2364485
export flank_kb=1000
export b1=6300000
export b2=6700000
export bracket=${b1}-${b2}

tabix ${INF}/METAL/gwas2vcf/TNFB.tsv.gz ${chr}:${bracket} > TNFB.tbx

module load python/2.7
awk -vchr=${chr} -vpos=${pos} -vd=$((${flank_kb}*1000)) 'BEGIN{print chr,pos-d,pos+d}' > st.tmp

INTERVAL
eQTLGen
SCALLOP
cd -

function hyprcoloc()
{
  join -a2 -e "NA" -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6,2.7 \
       <(sed '1d' ${INF}/MS/ieu-a-1025.assoc | \
         awk '{
                 chr=$2;pos=$4;A1=toupper($9);A2=toupper($10);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$3/$6,A1,A2,$8
              }' | \
         sort -k1,1 \
        ) \
       <(sed '1d' ${INF}/work/eQTLGen.lz | \
         awk '{
                 chr=$3;pos=$4;A1=toupper($5);A2=toupper($6);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$7,$3,$4,A1,A2,$2
              }' | \
         sort -k1,1 \
         ) | \
  join -16 \
       - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10/$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($4!=$10) $7=-$7
    if($4!=$13) $13=-$13
    print $1,$6,$8,$9,$4,$5,$3,$7,$13
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/work/LTBR.gassoc

  cut -d' ' -f1 ${INF}/work/LTBR.gassoc > ${INF}/work/LTBR.snpid
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/work/LTBR.snpid --r square --out ${INF}/work/LTBR
}

R --no-save -q <<END
  INF <- Sys.getenv("INF")
  library(gassocplot)
  d <- read.table(file.path(INF,"work","LTBR.gassoc"),col.names=c("snpid","marker","chr","pos","A1","A2","MS","LTBR","TNFB"))
  markers <- d[c("marker","chr","pos")]
  ld <- read.table(file.path(INF,"work","LTBR.ld"),col.names=with(d,marker),row.names=with(d,marker))
  z <- d[c("MS","LTBR","TNFB")]
  rownames(z) <- with(d,marker)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("MS","LTBR","TNFB"), ylab = "-log10(P)", legend=TRUE)
  pdf(file.path(INF,"plots","LTBR.pdf"),width=8,height=13)
  grid::grid.draw(sap)
  dev.off()
  stack_assoc_plot_save(sap, "LTBR.png", 5, width=8, height=13, dpi=300)
END

# ieu-a-1025.assoc
# p	chr	beta	position	n	se	id	rsid	ea	nea	eaf	trait
# eQTLGen.lz
# Pvalue	SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR	BonferroniP
# TNFB.lz
# MarkerName	Chromosome	Position	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Effect	StdErr	P	Direction	HetISq	HetChiSq	HetDf	logHetP	N

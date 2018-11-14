#!/bin/bash
. /etc/profile.d/modules.sh

# General notes, 14/11/18 JHZ
# 1. The overall design considers the fact that snpid (chr:pos_a1_a2) instead of rsid is used in the metal-analysis.
# 2. The snpid-rsid correspondence is obtained from snpstats_typed() and snpstats_imputed(), respectively.
# 3. PLINK clumping (clumped) provides corroborative result to GCTA -cojo (jma) used for PhenoScanner|cis/trans expliotation.
# 4. SNP-gene matchings are established by snp_gene() genomewide and genic_cis_trans() and Jimmy's R code for cis/trans classification.
#    Addtional notes:
#    - This follows https://github.com/jinghuazhao/PW-pipeline/blob/master/vegas2v2.sh
#    - bedtools 2.4.26 on cardio has no intersect command:
#    - module load bedtools/2.4.26
#    - We then compiled the latest bedtools release 2.27.1 to /scratch/jhz22/bin and gcc/4.8.1 is customarily called.
#    - The breakup of snpid leads to duplicate records in BED files so we employ uniq operation.
# 5. format_for_METAL() is actually copied from format.sh and CD6() is a simple exposition.

export BGEN_DIR=/scratch/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
export BGEN=$BGEN_DIR/interval_olink_subset_unrelated_4994_pihat_0.1875_autosomal_typed_only
function snpstats_typed()
{
  qctool -g $BGEN.bgen -s $BGEN.sample -snp-stats -osnp INTERVAL.snpstats
  awk 'NR==10' INTERVAL.snpstats | \
  awk '{gsub(/\t/, "\n",$0)};1'| \
  awk '{print "#" NR, $1}'
  cut -f1,2,4,5,6 INTERVAL.snpstats | \
  awk '(NR>10 && !/success/)' | \
  awk '
  {
     CHR=$1
     POS=$3
     a1=$4
     a2=$5
     if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
     else snpid="chr" CHR ":" POS "_" a1 "_" a2
     print snpid, $2
  }' | \
  sort -k1,1 > INTERVAL.snpid
}
# list of columns from operations above
#1 alternate_ids
#2 rsid
#3 chromosome
#4 position
#5 alleleA
#6 alleleB
#7 comment
#8 HW_exact_p_value
#9 HW_lrt_p_value
#10 alleleA_count
#11 alleleB_count
#12 alleleA_frequency
#13 alleleB_frequency
#14 minor_allele_frequency
#15 minor_allele
#16 major_allele
#17 info
#18 impute_info
#19 missing_proportion
#20 A
#21 B
#22 AA
#23 AB
#24 BB
#25 NULL
#26 total

export REF=/scratch/curated_genetic_data/reference_files/interval/
export TMPDIR=/scratch/jhz22/INF/work
function snpstats_imputed()
{
  cd work
  seq 22 | \
  parallel -j5 --env REF -C' ' '
    cut -f2-6 $REF/impute_{}_interval.snpstats | \
    awk "NR>1" | \
    awk "{\
      rsid=\$1;chr=\$2;pos=\$3;a1=\$4;a2=\$5; \
      gsub(/^0/,\"\",chr); \
      if (a1>a2) snpid=\"chr\" chr \":\" pos \"_\" a2 \"_\" a1; \
      else snpid=\"chr\" chr \":\" pos \"_\" a1 \"_\" a2; \
      print snpid, rsid \
    }" | \
    gzip -f > INTERVAL.snpid-{}.gz'
  cat INTERVAL.snpid-1.gz \
      INTERVAL.snpid-2.gz \
      INTERVAL.snpid-3.gz \
      INTERVAL.snpid-4.gz \
      INTERVAL.snpid-5.gz \
      INTERVAL.snpid-6.gz \
      INTERVAL.snpid-7.gz \
      INTERVAL.snpid-8.gz \
      INTERVAL.snpid-9.gz \
      INTERVAL.snpid-10.gz \
      INTERVAL.snpid-11.gz \
      INTERVAL.snpid-12.gz \
      INTERVAL.snpid-13.gz \
      INTERVAL.snpid-14.gz \
      INTERVAL.snpid-15.gz \
      INTERVAL.snpid-16.gz \
      INTERVAL.snpid-17.gz \
      INTERVAL.snpid-18.gz \
      INTERVAL.snpid-19.gz \
      INTERVAL.snpid-20.gz \
      INTERVAL.snpid-21.gz \
      INTERVAL.snpid-22.gz > INTERVAL.snpid.gz
  rm INTERVAL.snpid-*.gz
  cd -
}

function clumped_jma()
{
  cd work
  sed 's|work/INTERVAL.||g;s/.clumped://g' INTERVAL.clumped | \
  awk '{$1=$1;if(NR==1) $1="prot";if(NF>1) print}' > INTERVAL.clumped.dat

  sed 's|work/INTERVAL.||g;s/.ldr.cojo:/ /g' INTERVAL.ldr| \
  awk '{$1=$1; if(NR>1 && NF>1) print}' > INTERVAL.ldr.dat

  sed 's|work/INTERVAL.||g;s/.jma.cojo:/ /g' INTERVAL.jma | \
  awk '{$1=$1;if(NR==1) $1="prot";if(NF>1) print}' > INTERVAL.jma.dat
  awk 'NR>1' INTERVAL.jma.dat | \
  sort -k2,2n -k4,4n | \
  awk '{print "chr" $2 ":" $4}' | \
  uniq > INTERVAL.ps
# see also http://www.phenoscanner.medschl.cam.ac.uk
  module load phenoscanner/phenoscanner_v1.1
  phenoscanner -c All -l No -p 0.00001 -i INTERVAL.ps -o INTERVAL
  awk 'NR>1' INTERVAL.jma.dat | \
  sort -k3,3 | \
  cut -d' ' -f1,3 > INTERVAL.jma.snpid
  gunzip -c INTERVAL.snpid.gz | \
  sort -k1,1 | \
  join -11 -22 - INTERVAL.jma.snpid > INTERVAL.snpid_rsid
  cut -d' ' -f2 INTERVAL.snpid_rsid | \
  sort | \
  uniq > INTERVAL.rsid
  cd -
}

module load gcc/4.8.1
export INF=/scratch/jhz22/INF
function snp_gene()
# genomwide SNP-gene matchings
{
  cd work
  mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e 'select * from refGene' > refGene.txt
  cut -f3,7,8,13 refGene.txt | \
  awk '!index($1,"_")' | \
  uniq > refGene.bed
  cut -f3,7 $INF/doc/olink.inf.panel.annot.tsv  | \
  awk '(NR>1){
    gsub(/\"/,"",$0)
    if($1=="Q8NF90") $2="FGF5"
    if($2=="Q8WWJ7") $2="CD6"
    print $2
  }' | \
  sort > olink.gene
  (
    head -1 refGene.bed
    awk 'NR>1' refGene.bed | \
    sort -k4,4 | \
    join -11 -24 olink.gene - | \
    awk '{print $2,$3,$4,$1}'
  ) > refGene.olink
  wget -qO- https://www.cog-genomics.org/static/bin/plink/glist-hg19 > glist-hg19
  sort -k1,1n -k2,2n glist-hg19 | \
  awk '{if(NR==1) print "#chrom","start","end","gene";print "chr" $1,$2,$3,$4}' OFS="\t" > glist-hg19.bed
  (
    head -1 glist-hg19.bed
    awk 'NR>1' glist-hg19.bed | \
    sort -k4,4 | \
    join -11 -24 olink.gene - | \
    awk '{print $2,$3,$4,$1}'
  ) > glist-hg19.olink
  awk -vOFS="\t" '{
    snpid=$1
    rsid=$2
    split(snpid,a,":")
    chr=a[1]
    split(a[2],b,"_")
    pos=b[1]
    prot=$3
    if(NR==1) print "#chrom","start","end","rsid","prot"
    print chr,pos-1,pos,rsid,prot
  }' INTERVAL.snpid_rsid > INTERVAL.bed
  (
    echo -e "#chrom\tstart\tend\t\tgene\tprot";
    sort -k2,2 $INF/inf1.list > inf1.tmp
    awk '{
      FS=OFS="\t"
      gsub(/\"/,"",$0)
      if($3=="Q8NF90") $7="FGF5"
      if($3=="Q8WWJ7") $7="CD6"
      print "chr" $8,$9,$10,$3,$7
    }' $INF/doc/olink.inf.panel.annot.tsv | \
    sort -k4,4 | \
    join -14 -22 - inf1.tmp | \
    awk -vOFS="\t" '{print $2,$3,$4,$5,$6}'
  ) > olink.bed
  awk -vOFS="\t" -vM=1000000 '{
    chrom=$1
    cdsStart=$2
    cdsEnd=$3
    gene=$4
    prot=$5
    start=cdsStart-M
    if (start<0) start=0
    end=cdsEnd+M
    if(NR==1) print "#chrom", "start", "end", "gene", "prot";
    else print chrom, start, end, gene, prot
  }' olink.bed > olink.cis
  bedtools intersect -a INTERVAL.bed -b refGene.bed -loj > INTERVAL.refGene
  bedtools intersect -a INTERVAL.bed -b glist-hg19.bed -loj > INTERVAL.glist-hg19
  bedtools intersect -a INTERVAL.bed -b olink.bed -loj > INTERVAL.olink
  bedtools intersect -a INTERVAL.bed -b olink.cis -loj > INTERVAL.genic_cis_trans
  cd -
}
#  head -1 $INF/doc/olink.inf.panel.annot.tsv | \
#  awk '{gsub(/\t/, "\n",$0)};1'| \
#  awk '{print "#" NR, $1}'
#1 "target"
#2 "target.short"
#3 "uniprot"
#4 "panel"
#5 "prot.on.multiple.panel"
#6 "panels.with.prot"
#7 "hgnc_symbol"
#8 "chromosome_name"
#9 "start_position"
#10 "end_position"
#11 "olink.id"
#12 "alternate.uniprot"

R --no-save -q < $INF/doc/cis.vs.trans.classification.R

function genic_cis_trans()
# title,genic,cis,trans
{
  awk '$NF!="."' INTERVAL.olink > INTERVAL.genic
  cut -f4 INTERVAL.genic > INTERVAL.rsid_genic
  grep -v -w -f INTERVAL.rsid_genic INTERVAL.bed > INTERVAL.tmp
  bedtools intersect -a INTERVAL.tmp -b olink.cis -loj > INTERVAL.cis_trans
  awk '$5==$10' INTERVAL.cis_trans > INTERVAL.cis
  cut -f4 INTERVAL.cis > INTERVAL.rsid_cis
  awk '$5!=$10&&$10=="."' INTERVAL.cis_trans > INTERVAL.trans
  cut -f4 INTERVAL.trans > INTERVAL.rsid_trans
  wc -l INTERVAL.rsid_cis INTERVAL.rsid_genic INTERVAL.rsid_trans
  sort -k5,5 INTERVAL.trans | \
  join -t$'\t' -11 -25 $INF/inf1.list - | \
  sort -k2,2 > olink.tmp
  cut -f3,7 $INF/doc/olink.inf.panel.annot.tsv | \
  awk -vOFS="\t" 'NR>1{
    gsub(/\"/,"",$0);
    if($1=="Q8NF90") $2="FGF5"
    if($1=="Q8WWJ7") $2="CD6"
    print $1,$2
  }' | \
  sort -k1,1 | \
  join -t$'\t' -12 -21 olink.tmp - | \
  awk -vOFS="\t" '{print $3,$4,$5,$6,$2,$8,$9,$10,$12,$11}' > INTERVAL.tmp
  (
    echo -e "chr\tstart\tend\trsid\tprot\tchrom\tcisstart\tcisend\tgene\ttarget\tstatus"
    awk -vOFS="\t" '{print $0,"genic"}' INTERVAL.genic
    awk -vOFS="\t" '{print $0,"cis"}' INTERVAL.cis
    awk -vOFS="\t" '{print $0,"trans"}' INTERVAL.tmp
  ) > INTERVAL.prot
  R --no-save -q <<END
    prot <- read.delim("INTERVAL.prot",as.is=TRUE,header=TRUE)
    genic_cis_trans <- with(prot,as.matrix(table(gene,status)))
    genic_cis_trans
    sum(genic_cis_trans)
END
}

genic_cis_trans

export INTERVAL=/scratch/jp549/olink-merged-output
function format_for_METAL()
{
  ls $INTERVAL/*gz | \
  grep inf1 | \
  xargs -l basename | \
  sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
  cut -d'_' --output-delimiter=' ' -f1,4 | \
  parallel -j2 --env INTERVAL -C' ' '
    gunzip -c $INTERVAL/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
    awk -f files/INTERVAL.awk | \
    awk -f files/order.awk | \
    gzip -f > sumstats/INTERVAL.{1}.gz
  '
}

function CD6()
# SUMSTATS for depict
{
  gunzip -c $INTERVAL/INTERVAL_inf1_CD6___Q8WWJ7_chr_merged.gz | \
  awk -vOFS="\t" '(NR>1){
       SNP=$2;
       chr=$3; sub(/^0/,"",chr);
       pos=$4;
       a1=$5;
       a2=$6;
       N=$18;
       EAF=(0.5*$15+$16)/N;
       beta=$24;
       se=$25;
       p=$22;
       if(SNP!="." && p<=0.1) print SNP,a2,a1,EAF,beta,se,p,N,chr,pos;
  }' | \
  sort -k9,9n -k10,10n > CD6
}

export PHEN=/scratch/curated_genetic_data/phenotypes/interval/high_dimensional_data/Olink_proteomics_inf/gwasqc/olink_qcgwas_inf.csv
export IMPUTED=/scratch/curated_genetic_data/interval/imputed
function reference()
{
  cut -d"," -f1 $PHEN | \
  awk 'NR>1 {OFS="\t";print $1,$1}' > INTERVAL.id
  seq 22 | \
  parallel -j3 --env IMPUTED -C' ' '
    plink --bgen $IMPUTED/impute_{}_interval.bgen \
          --sample $IMPUTED/interval.samples \
          --keep INTERVAL.id \
          --make-bed \
          --out INTERVAL-{} \
          --threads 2
  '
}

export SCRIPT=/scratch/jp549/analyses/interval_subset_olink/inf1/r2/outlier_in/pcs1_3
export BS=/scratch/jp549/apps/bram-scripts

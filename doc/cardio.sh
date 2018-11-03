#!/bin/bash
. /etc/profile.d/modules.sh

# 3-11-2018 JHZ

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

function process_clump_cojo()
{
  cd work
  sed 's|work/INTERVAL.||g;s/.clumped://g' INTERVAL.clumped | \
  awk '{$1=$1;if(NR==1) $1="prot";if(NF>1) print}' > INTERVAL.clumped.dat

  sed 's|work/INTERVAL.||g;s/.jma.cojo:/ /g' INTERVAL.jma | \
  awk '{$1=$1;if(NR==1) $1="prot";if(NF>1) print}' > INTERVAL.jma.dat
  awk 'NR>1' INTERVAL.jma.dat | \
  sort -k2,2n -k4,4n | \
  awk '{print "chr" $2 ":" $4}' | \
  uniq > INTERVAL.ps

  sed 's|work/INTERVAL.||g;s/.ldr.cojo:/ /g' INTERVAL.ldr| \
  awk '{$1=$1; if(NR>1 && NF>1) print}' > INTERVAL.ldr.dat
  cd -
}

function signals()
{
  cd work
  awk 'NR>1' INTERVAL.jma.dat | \
  sort -k3,3 | \
  cut -d' ' -f1,3 > INTERVAL.jma.snpid
  gunzip -c INTERVAL.snpid.gz | \
  sort -k1,1 | \
  join -11 -22 - INTERVAL.jma.snpid > INTERVAL.snpid_rsid
  cut -d' ' -f2 INTERVAL.snpid_rsid > INTERVAL.rsid
  cd -
}

function ps()
# see also http://www.phenoscanner.medschl.cam.ac.uk
{
  module load phenoscanner/phenoscanner_v1.1
  cd work
  phenoscanner -c All -l No -p 0.00001 -i INTERVAL.ps -o INTERVAL
  cd -
}

export INTERVAL=/scratch/jp549/olink-merged-output
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

function snp_gene()
{
  cd work
  mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e 'select * from refGene' > refGene.txt
  # https://github.com/jinghuazhao/PW-pipeline/blob/master/vegas2v2.sh 
  wget https://www.cog-genomics.org/static/bin/plink/glist-hg19
  sort -k1,1n -k2,2n glist-hg19 | \
  awk '{if(NR==1) print "#chrom","start","end","gene";print "chr" $1,$2,$3,$4}' OFS="\t" > glist-hg19.bed
  awk '{
    snpid=$1
    rsid=$2
    split(snpid,a,":")
    chr=a[1]
    split(a[2],b,"_")
    pos=b[1]
    if(NR==1) print "#chrom","Start","End","rsid"
    print chr,pos-1,pos,rsid
  }' INTERVAL.snpid_rsid > INTERVAL.bed
# The following module is available on cardio but it does not contain the command.
# module load bedtools/2.4.26
  intersectBed -a INTERVAL.bed -b glist-hg19.bed -loj > INTERVAL.bedtools
  cd -
}

export PHEN=/scratch/curated_genetic_data/phenotypes/interval/high_dimensional_data/Olink_proteomics_inf/gwasqc/olink_qcgwas_inf.csv
export SCRIPT=/scratch/jp549/analyses/interval_subset_olink/inf1/r2/outlier_in/pcs1_3
export BS=/scratch/jp549/apps/bram-scripts
function NOTE_for_MASS_Processing()
{
  ls $INTERVAL/*gz | \
  grep inf1 | \
  xargs -l basename | \
  sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
  cut -d'_' --output-delimiter=' ' -f1,4 | \
  grep -w CD6 | \
  parallel -j2 --env INTERVAL -C' ' '
    gunzip -c $INTERVAL/INTERVAL_inf1_{1}___{2}_chr_merged.gz
  '
}

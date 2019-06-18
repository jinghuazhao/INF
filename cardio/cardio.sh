#!/bin/bash
. /etc/profile.d/modules.sh

# General notes, 29/5/19 JHZ
# 1. The overall design considers the fact that snpid (chr:pos_a1_a2) instead of rsid is used in the metal-analysis.
# 2. The snpid-rsid correspondence is obtained from snpstats_typed() and snpstats_imputed(), respectively.
# 3. PLINK clumping (clumped) provides corroborative result to GCTA -cojo (jma) used for PhenoScanner|cis/trans expliotation.
# 4. SNP-gene matchings are established by snp_gene() genomewide and Jimmy's R code for cis/trans classification.
#    Addtional notes:
#    - This follows https://github.com/jinghuazhao/PW-pipeline/blob/master/vegas2v2.sh
#    - bedtools 2.4.26 on cardio has no intersect command:
#    - module load bedtools/2.4.26
#    - We then compiled the latest bedtools release 2.27.1 to /scratch/jhz22/bin and gcc/4.8.1 is customarily called.
#    - The breakup of snpid leads to duplicate records in BED files so we employ uniq operation.
# 5. format_for_METAL() is actually copied from format.sh.

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

export TMPDIR=/scratch/jhz22/tmp

(
  seq 22 | \
  parallel -j1 --env REF -C' ' '
    cat $REF/impute_{}_interval.snpstats | \
    awk -vOFS="\t" "{
      if(NR==1) print \"SNPID\",\"N\",\"MAF\",\"HWE\",\"info\";
      else {
        CHR=\$3
        sub(/^0/,\"\",CHR)
        POS=\$4
        a1=\$5
        a2=\$6
        N=\$9+\$10+\$11
        MAF=\$15
        HWE=\$16
        info=\$19
        if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
        else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
        print snpid, N, MAF, HWE, info
      }
    }"
  '
) | \
awk '(NR==1||$1!="SNPID")' | \
sort -k1,1 | \
gzip -f > INTERVAL.snpstats.gz

# /scratch/curated_genetic_data/reference_files/interval/
# gunzip -c impute_1_interval.snpstats.gz | \
# head -1 | \
# awk '{gsub(/\t/, "\n",$0)};1'| awk '{print "#" NR, $1}'

#1 SNPID
#2 RSID
#3 chromosome
#4 position
#5 A_allele
#6 B_allele
#7 minor_allele
#8 major_allele
#9 AA
#10 AB
#11 BB
#12 AA_calls
#13 AB_calls
#14 BB_calls
#15 MAF
#16 HWE
#17 missing
#18 missing_calls
#19 information

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

R --no-save -q < $INF/cardio/cis.vs.trans.classification.R

export INTERVAL=/scratch/jp549/olink-merged-output
function format_for_METAL()
{
  ls $INTERVAL/*gz | \
  grep inf1 | \
  xargs -l basename | \
  sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g' | \
  cut -d'_' --output-delimiter=' ' -f1,4 | \
  parallel -j6 --env INTERVAL -C' ' '
    gunzip -c $INTERVAL/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
    awk -f tryggve/INTERVAL.awk | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/INTERVAL.{1}.gz
  '
}

echo "--> clumping"

export rt=$HOME/INF/METAL
sbatch --wait cardio/clump.sb
(
  grep CHR $rt/*.clumped | \
  head -1
  grep -v CHR $rt/*.clumped
) | \
sed 's|'"$rt"'/||g;s/.clumped://g' | \
awk '(NF>1){$3="";print}' | \
awk '{$1=$1;if(NR==1)$1="prot";print}' > INF1.clumped
R --no-save -q <<END
  require(gap)
  clumped <- read.table("INF1.clumped",as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","BP","SNP","prot")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","uniprot")
  cistrans <- cis.vs.trans.classification(hits)
  sink("INF1.clumped.out")
  with(cistrans,table)
  sink()
  sum(with(cistrans,table))
  pdf("INF1.circlize.pdf")
  circos.cis.vs.trans.plot(hits="INF1.clumped")
  dev.off()
END

echo "--> AILD collection"

function aild_clump()
{
  if [ -f INF1.clumped ]; then rm INF1.clumped; fi
  (
    cat *.clumped | \
    head -1 | awk '{print "prot", $0}'
    for i in $(ls *clumped); do awk -vi=$i 'NR>1{split(i,a,"-");print a[1],$0}' $i; done
  ) | \
  sed 's/.clumped://g' | \
  awk '(NF>1){$3="";print}' | \
  awk '{$1=$1;if(NR==1)$1="prot";print}' > INF1.clumped
}

function aild_cojo()
{
  if [ -f INF1.jma ]; then rm INF1.jma; fi
  (
    cat *.jma.cojo | \
    head -1 | awk '{print "prot", $0}'
    for i in $(ls *jma.cojo); do awk -vi=$i 'NR>1{split(i,a,"-");print a[1],$0}' $i; done
  ) | \
  sed 's/.jma.cojo://g' | \
  awk '{$1=$1;if(NR==1)$1="prot";print}'| sed 's/ /\t/g' > INF1.jma
}

echo "--> finemapping"

module load gcc/5.2.0

(
  echo "chr start end snpid pos r prot"
  (
    echo -e "chrom start end SNPID pos prot" 
    awk 'NR>1 {print "chr" $2, $4-1, $4, $3, $4, $1}' INF1.clumped 
  ) | \
  sed 's/ /\t/g' | \
  bedtools intersect -a /scratch/jhz22/FM-pipeline/1KG/EUR.bed -b - -loj | \
  awk '$5!="." {
    gsub(/chr/,"",$1);
    gsub(/region/,"",$4);
    $5="";$6="";$7="";
    print $1,$2,$3,$8,$9,$4,$10}'
) > st.bed

awk 'NR>1 {print $1,$2,$3}' INF1.clumped | \
parallel -j3 -C' ' '
  gunzip -c METAL/{1}-1.tbl.gz | \
  awk -vchr={2} "chr==\$1{print \$3,toupper(\$4),toupper(\$5),\$6,\$10,\$11,\$12,\$14,\$1,\$2}" > METAL/{1}-{3}
'

export wd=${PWD}
export rt=$wd/METAL/
source $wd/cardio/fm.ini
export FM_location=$HOME/FM-pipeline

awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' \
         --env wd \
         --env rt \
         --env FM_location \
         --env GEN_location  \
         --env CAVIAR \
         --env CAVIARBF \
         --env FM_summary \
         --env fgwas \
         --env finemap \
         --env GCTA \
         --env JAM \
         --env LD_MAGIC \
         --env LD_PLINK \
         --env LocusZoom \
          '$wd/cardio/fm.subs {1} {2} {3} {4} {5} {6} {7}'

## finemapping cojo results

export rt=snps/cojo
(
  echo "chr start end snpid pos r prot"
  (
    echo -e "chrom start end SNPID pos prot"
    awk 'NR>1{print "chr" $2, $4-1, $4, $3, $4, $1}' $rt/INF1.jma 
  ) | \
  sed 's/ /\t/g' | \
  bedtools intersect -a /scratch/jhz22/FM-pipeline/1KG/EUR.bed -b - -loj | \
  awk '$5!="." {
    gsub(/chr/,"",$1);
    gsub(/region/,"",$4);
    $5="";$6="";$7="";
    print $1,$2,$3,$8,$9,$4,$10}'
) > st.bed

awk 'NR>1 {print $1,$2,$3}' $rt/INF1.jma | \
parallel -j3 -C' ' '
  gunzip -c METAL/{1}-1.tbl.gz | \
  awk -vchr={2} "chr==\$1{print \$3,toupper(\$4),toupper(\$5),\$6,\$10,\$11,\$12,\$14,\$1,\$2}" > work/{1}-{3}
'

export wd=${PWD}
export rt=$wd/work/
source $wd/cardio/fm.ini
export FM_location=$HOME/FM-pipeline

awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' \
         --env wd \
         --env rt \
         --env FM_location \
         --env GEN_location  \
         --env CAVIAR \
         --env CAVIARBF \
         --env FM_summary \
         --env fgwas \
         --env finemap \
         --env GCTA \
         --env JAM \
         --env LD_MAGIC \
         --env LD_PLINK \
         --env LocusZoom \
          '$wd/cardio/fm.subs {1} {2} {3} {4} {5} {6} {7}'

ls /home/jhz22/INF/sumstats/INTERVAL/*gz | \
sed 's|/home/jhz22/INF/sumstats/INTERVAL/INTERVAL.||g;s/.gz//g' > manhattan.list

sbatch manhattan.sb

# distcojo
(
  cat *jma.cojo | \
  head -1 | \
  awk -vOFS="\t" '{print "prot","MarkerName",$0}'
  ls *jma.cojo | \
  sed 's/-/ /g;s/.jma.cojo//g' | \
  parallel -j1 -C' ' 'awk -vOFS="\t" -vprot={1} -vMarkerName={2} "NR > 1{print prot,MarkerName,\$0}" {1}-{2}.jma.cojo'
) > INF1.jma

# 1-11-2018 JHZ

export INTERVAL=/scratch/jp549/olink-merged-output

function CD6()
{
  gunzip -c $INTERVAL/INTERVAL_inf1_CD6___Q8WWJ7_chr_merged.gz | \
  awk -vOFS="\t" '(NR>1){
       SNP=$2;
       chr=$3; sub(/0/,"",chr);
       pos=$4;
       a1=$5;
       a2=$6;
       N=$18;
       EAF=(0.5*$15+$16)/N;
       beta=$24;
       se=$25;
       p=$22;
       if(p<=1e-5) print SNP,a2,a1,EAF,beta,se,p,N,chr,pos;
  }' | \
  sort -k9,9n -k10,10n > CD6
}

export BGEN_DIR=/scratch/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_typed_only
export BGEN=$BGEN_DIR/interval_olink_subset_unrelated_4994_pihat_0.1875_autosomal_typed_only

function snpstats()
{
  qctool -g $BGEN.bgen -s $BGEN.sample -snp-stats -osnp INTERVAL.snpstats

  awk 'NR==10' INTERVAL.snpstats | \
  awk '{gsub(/\t/, "\n",$0)};1'| \
  awk '{print "#" NR, $1}'
}

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

sed 's|work/INTERVAL.||g;s/.clumped://g' work/INTERVAL.clumped| \
awk '{$1=$1;if(NR==1) $1="prot";if(NF>1) print}' > work/INTERVAL.clump.dat

function signals()
{
  head -1 INTERVAL.snpstats
  awk '
  {
     OFS="\t"
     if (NR>1)
     {
       CHR=$3
       POS=$4
       a1=$5
       a2=$6
       if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
       else snpid="chr" CHR ":" POS "_" a1 "_" a2
       $1=snpid
     }
     print snpid, $2
  }' INTERVAL.snpstats | \
  sort | \
  join - INTERVAL.clumped | \
  cut -d' ' -f2 > INTERVAL.rsid
}

export SCRIPT=/scratch/jp549/analyses/interval_subset_olink/inf1/r2/outlier_in/pcs1_3
export BS=/scratch/jp549/apps/bram-scripts
export INF=/scratch/curated_genetic_data/phenotypes/interval/high_dimensional_data/Olink_proteomics_inf/gwasqc/olink_qcgwas_inf.csv
export REF=/scratch/curated_genetic_data/reference_files/interval/

function NOTE()
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

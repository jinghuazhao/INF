# 1-4-2020 JHZ

export threads=6

function ARISTOTLE()
{
  export rt=/data/niceri/ARISTOTLE_INFL/infl_rnorm_lod_for_meta
  ls ${rt} | \
  sed 's/ARISTOTLE.infl.//g;s/.20200120.chr/ /g;s/.txt.gz//g' | \
  cut -d' ' -f1 | \
  uniq | \
  parallel -j${threads} -C' ' --env rt '
  (
     gunzip -c ${rt}/ARISTOTLE.infl.{}.20200120.chr1.txt.gz | \
     head -1;
     for i in $(seq 22)
     do
       gunzip -c ${rt}/ARISTOTLE.infl.{}.20200120.chr${i}.txt.gz | \
       awk -f tryggve/ARISTOTLE.awk | \
       awk "a[\$1]++==0" | \
       sort -k2,2n -k3,3n -t$''\t''
     done
  ) | \
  gzip -f > ARISTOTLE.{}.txt.gz
 '
}

function BioFinder()
{
  grep -v -w TNF sumstats/BioFinder.list | \
  parallel -j$threads -C' ' '
     awk -f tryggve/BioFinder.awk /data/andmala/biofinder_inf/rsannot_runGwas_plasmaImp.{1}_zre_INFI.glm.linear | \
     awk -f tryggve/order.awk | \
     gzip -f > sumstats/BioFinder/BioFinder.{3}.gz'
# version with complete data
  gunzip -c /data/jinhua/data/BioFinder/rsannot_runGwas_plasmaImp.TNF_zre_INFI.glm.linear.gz | \
  awk -f tryggve/BioFinder.awk | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/BioFinder/BioFinder.TNF.gz
}

function EDCUT()
{
# SNPID has prefix esv for non-rsids
  cat sumstats/EGCUT.list | \
  parallel -j$threads -C' ' '
    gunzip -c /data/anekal/EGCUT_INF/EGCUT_autosomal_{1}_inf_280918.txt.gz | \
    awk "{if(NR>1&&(index(\$1,\"esv\")||index(\$1,\"ss\"))) \$1=\"chr\" \$2 \":\" \$3;if(\$13>0.4) print}" | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/EGCUT/EGCUT.{2}.gz'
}

function INTERVAL()
{
  cat sumstats/INTERVAL.list | \
  parallel -j$threads -C' ' '
     /usr/bin/gunzip -c /data/jampet/upload-20170920/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
     awk -f tryggve/INTERVAL.awk | \
     awk -f tryggve/order.awk | \
     gzip -f > sumstats/INTERVAL/INTERVAL.{1}.gz'
}

function LifeLinesDeep()
{
# SNPID has no "chr" prefix for non-rsids
  cat sumstats/LifeLinesDeep.list | \
  sed 's/_/\t/g' | \
  cut -f2 | \
  sort | \
  join -11 -23  - inf1_gene | \
  parallel -j$threads -C' ' 'gunzip -c /data/darzhe/LifeLinesDeep.cistranspQTLs.20171220.txt.gz | \
  awk -vprotein={1} -vFS="\t" -vOFS="\t" "(NR==1||index(\$1,protein))" | \
  cut -f2-14 | \
  sort -k2,2n -k3,3n | \
  awk "{if (NR>1&&substr(\$1,1,2)!=\"rs\") \$1=\"chr\" \$2 \":\" \$3; print}" | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/LifeLinesDeep/LifeLinesDeep.{1}.gz'
}

function KORA()
{
# 91 proteins without BDNF P23560 BDNF
  cat sumstats/KORA.list | \
  parallel -j$threads -C' ' '
     gunzip -c KORA/{1}.gz | \
     awk -f tryggve/KORA.awk | \
     awk -f tryggve/order.awk | \
     gzip -f > sumstats/KORA/KORA.{3}.gz'
}

function MadCam()
{
  cat sumstats/MadCam.list | \
  parallel -j$threads -C' ' '
     cut -f1-14 /data/andmala/madcam/MadCAM.{1}.{2}.txt | \
     sed 's/CODE_ALLELE_FQ/CODE_ALL_FQ/g' | \
     awk "\$13>0.3" | \
     awk -f tryggve/order.awk | \
     gzip -f > sumstats/MadCam/MadCam.{3}.gz'
}

function NSPHS()
{
  export NSPHS=/data/stefane/NSPHS_INF
  cat sumstats/NSPHS.list | \
  parallel -j$threads --env NSPHS -C' ' '
    gunzip -c $NSPHS/NSPHS_inf1_{3}_{1}.txt.gz | \
    awk -f tryggve/NSPHS.awk | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/NSPHS/NSPHS.{2}.gz'
}

function PIVUS_ULSAM()
{
# SNPID has :I/D suffix and VG prefix
  ls /data/stefang/pivus_ulsam/pivus* | \
  xargs -l -x basename | \
  sed 's/pivus.all.//g;s/.20161128.txt.gz//g' | \
  sort | \
  join - inf1_gene | \
  cut -d ' ' -f1,3 | \
  parallel -j$threads -C' ' '
    gunzip -c /data/stefang/pivus_ulsam/pivus.all.{1}.20161128.txt.gz | \
    awk "{if(NR>1&&substr(\$1,1,2)!=\"rs\") \$1=\"chr\" \$2 \":\" \$3;print}" | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/PIVUS/PIVUS.{2}.gz'

  ls /data/stefang/pivus_ulsam/ulsam* | \
  xargs -l -x basename | sed 's/ulsam.all.//g;s/.20161128.txt.gz//g' | \
  sort | \
  join - inf1_gene | \
  cut -d ' ' -f1,3 | \
  parallel -j$threads -C' ' '
    gunzip -c /data/stefang/pivus_ulsam/ulsam.all.{1}.20161128.txt.gz | \
    awk "{if(NR>1&&substr(\$1,1,2)!=\"rs\") \$1=\"chr\" \$2 \":\" \$3;print}" | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/ULSAM/ULSAM.{2}.gz'
}

function ORCADES()
{
  cat sumstats/ORCADES.list | \
  parallel -j$threads -C' ' '
    gunzip -c /data/erimac/ORCADES/ORCADES.INF1.{1}_rank.tsv.gz | \
    awk "NR==1||\$13>0.4" | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/ORCADES/ORCADES.{2}.gz'
}

function VIS()
{
  cat sumstats/VIS.list | \
  parallel -j$threads -C' ' '
    gunzip -c /data/erimac/VIS/VIS.INF1.{1}_rank.tsv.gz | \
    awk "NR==1||\$13>0.4" | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/VIS/VIS.{2}.gz'
}

function RECOMBINE()
{
  export rt=/data/jinhua/data/RECOMBINE/RECOMBINE_pQTLs__meta_scallop
  export rt=/data/jinhua/data/RECOMBINE/RECOMBINE_INF1_pQTLs_updated_13thMarch_19
  sort -k2,2 inf1.list > inf1.tmp
  sort -k3,3 sumstats/RECOMBINE.list | \
  join -13 -22 - inf1.tmp | \
  parallel -j$threads --env rt -C' ' '
    (
      for chr in `seq 22`
      do
        if [ $chr -eq 1 ]; then
           gunzip -c $rt/{2}_{3}___{1}_chr${chr}_RECOMBINE.txt.gz | \
           awk "NR == 1"
        fi
        gunzip -c $rt/{2}_{3}___{1}_chr${chr}_RECOMBINE.txt.gz | \
        awk "NR > 1" 
      done
    ) | \
    awk "NR==1||(\$16>0.3){sub(/EFFECT_ALL_FQ/,\"CODE_ALL_FQ\",\$11);print}" | \
    sed "s/ /\t/g" | \
    cut -f3-7,9-17 | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/RECOMBINE/RECOMBINE.{4}.gz'
}

function STABILITY()
{
  export STABILITY=/data/niceri/Stability_INF1
  sort -k3,3 sumstats/STABILITY.list | \
  join -13 -21 - work/STABILITY.N | \
  parallel -j$threads --env STABILITY -C' ' '
    (
      for chr in `seq 22`; do zgrep -v -w SNPID $STABILITY/STABILITY_{2}_{3}_chr${chr}.txt.gz; done
    ) | \
    awk -vOFS="\t" -vN={4} -f tryggve/STABILITY.awk | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/STABILITY/STABILITY.{1}.gz'
}

function STANLEY_lah1()
{
  export STANLEY_lah1=/data/andmala/STANLEY_20180911
  export N=344
  cut -d' ' -f1-3 sumstats/STANLEY.list | \
  parallel -j$threads --env STANLEY_lah1 --env N -C' ' '
  (
    for chr in `seq 22`; do gunzip -c $STANLEY_lah1/STANLEY_lah1_inf_chr${chr}_pheno{1}.txt.assoc.dosage.gz; done
  ) | \
  awk "NR==1||\$2!=SNP" | \
  awk -vN=$N -f tryggve/STANLEY.awk | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/STANLEY/STANLEY_lah1.{2}.gz'
}

function STANLEY_swe6()
{
  export STANLEY_swe6=/data/andmala/STANLEY_20180911//swe6_inf
  export N=300
  cut -d' ' -f1-3 sumstats/STANLEY.list | \
  parallel -j$threads --env STANLEY_swe6 --env N -C' ' '
  (
    for chr in `seq 22`; do gunzip -c $STANLEY_swe6/STANLEY_swe6_inf_chr${chr}_pheno{1}.txt.assoc.dosage.gz; done
  ) | \
  awk "NR==1||\$2!=SNP" | \
  awk -vN=$N -f tryggve/STANLEY.awk | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/STANLEY/STANLEY_swe6.{2}.gz'
}

function checklines()
# to check for number of SNPs for all proteins in a particular study
# usage checklines KORA 3
# where the protein names are the 3rd column in KORA.list
{
  export study=$1
  export col=$2
  parallel --env study --env col -C' ' '
    echo $study-{}
    gunzip -c sumstats/$study/${study}.{}.gz | \
    wc -l' ::: $(cut -d" " -f$col sumstats/$study.list)
}

function turbo()
{
  parallel -j5 -C' ' '
    export s={1}
    export p={2}
    export g=$(grep -w {2} inf1.gene | cut -d" " -f3)
    zcat sumstats/${s}/${s}.${p}.gz | \
    awk "NR>1&&!/CHR/{print \$2,\$3,\$11}" | \
    gzip -f > ${s}.${p}.gz
    (
      zcat glist.gz | \
      head -1
      zgrep -w $g glist.gz
    ) | \
    gzip -f > ${s}.${p}.glist.gz
  # Manhattan
  R --slave --vanilla --args \
    input_data_path=${s}.${p}.gz \
    output_data_rootname=${s}.${p}.manhattan \
    custom_peak_annotation_file_path=${s}.${p}.glist.gz \
    reference_file_path=cardio/turboman_hg19_reference_data.rda \
    pvalue_sign=5e-10 \
    plot_title="Manhattan plot" < cardio/turboman.r
  # QQ
  R --slave --vanilla --args \
    input_data_path=${s}.${p}.gz \
    output_data_rootname=${s}.${p}.qq \
    plot_title="Q-Q plot" < cardio/turboqq.r
    rm ${s}.${p}.gz ${s}.${p}.glist.gz
  ' ::: INTERVAL BioFinder EGCUT MadCam KORA NSPHS ORCADES RECOMBINE STABILITY STANLEY VIS ::: $(cut -d' ' -f1 prot.list)
}

$1

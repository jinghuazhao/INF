# 19-2-2019 JHZ

# --- INF list of proteins and file list ---

function INF()
# INF list of proteins
{
  ## original list
  grep inf1 doc/olink.prot.list.txt | \
  sed 's/inf1_//g;s/___/\t/g' | \
  sort -k1,1 > inf1.list
  grep inf1 doc/olink.prot.list.txt | \
  sed 's/inf1_//g;s/___/\t/g' | \
  sort -k2,2 > prot.list
  awk -vOFS="\t" '{l=$1;gsub(/\./,"_",$1);print $1,$2,l}' inf1.list > inf1_gene
  ## adding aliases
  (
    grep inf1 doc/olink.prot.list.txt | \
    sed 's/inf1_//g;s/___/\t/g'
    echo -e "CD6\tP30203\nFGF.5\tP12034"
  ) > inf1.tmp
  sort -k1,1 inf1.tmp > inf1.list
  cut -f3,7 doc/olink.inf.panel.annot.tsv | \
  sed 's/\"//g' | \
  sort -k1,1 | \
  join -12 -21 prot.list - | \
  awk '{
       if($1=="Q8NF90") {$1="P12034";$3="TGF5"}
       if($1=="Q8WWJ7") {$1="P30203";$3="CD6"}
       print
  }' | \
  sort -k3,3 > inf1.gene
}

# --- studies ----

BioFinder()
{
  awk '
  {
    if ($3=="IL8") $3="CXCL8"
    if ($3=="TGF5") $3="FGF5"
  };1' inf1.gene | \
  sort -k3,3 > inf1.tmp
  ls /data/andmala/biofinder_inf/ | \
  grep -v list | \
  sed 's/rsannot_runGwas_plasmaImp.//g;s/_zre_INFI.glm.linear\*//g' | \
  sort -k1,1 | \
  join -11 -23 - inf1.tmp > sumstats/BioFinder.list
}

function EGCUT_INF() {
  # EGCUT_INF by autosomal, female, male
  ls /data/anekal/EGCUT_INF/ | \
  grep inf > sumstats/EGCUT_INF.list
}

function INTERVAL() {
  # INTERVAL, 92 lines
  ls /data/jampet/upload-20170920/ | \
  grep inf > sumstats/INTERVAL.list
}

function KORA() {
  awk '
  {
     if ($3=="EIF4EBP1") $3="4EBP1"
     if ($3=="NGF") $3="BetaNGF"
     if ($3=="S100A12") $3="ENRAGE"
     if ($3=="TGF5") $3="FGF5"
     if ($3=="FLT3LG") $3="FLT3L"
     if ($3=="IFNG") $3="IFNg"
     if ($3=="IL1A") $3="IL1a"
     if ($3=="TGFB1") $3="LAP_TGFb1"
     if ($3=="CCL2") $3="MCP1"
     if ($3=="CCL8") $3="MCP2"
     if ($3=="CCL7") $3="MCP3"
     if ($3=="CCL13") $3="MCP4"
     if ($3=="TNFRSF11B") $3="OPG"
     if ($3=="CD274") $3="PDL1"
     if ($3=="KITLG") $3="SCF"
     if ($3=="SULT1A1") $3="ST1A1"
     if ($3=="TGFA") $3="TGFa"
     if ($3=="LTA") $3="TNFB"
     if ($3=="TNFSF10") $3="TRAIL"
     if ($3=="TNFSF11") $3="TRANCE"
     if ($3=="TNFSF12") $3="TWEAK"
     if ($3=="PLAU") $3="uPA"
  };1' inf1.gene | \
  sort -k3,3 > inf1.tmp
  ls KORA/snptest.*.out.gz | \
  awk '{gsub(/KORA\/snptest.|.out.gz/,"");print}' | \
  sort -k1,1 | \
  join -11 -23 - inf1.tmp > sumstats/KORA.list
}

function LifeLines() {
  # LifeLines
  gunzip -c /data/darzhe/LifeLinesDeep.cistranspQTLs.20171220.txt.gz | \
  cut -f1 | \
  awk '(NR>1)' | \
  sort | \
  uniq > sumstats/LifeLinesDeep.list
}

function MadCam()
{
  ls /data/andmala/madcam | \
  awk '{split($1,a,".");print a[2],a[3]}' | \
  sort -k1,1 | \
  join -11 -22 - prot.list > sumstats/MadCam.list
}

function NSPHS_INF() {
  # NSPHS_INF
  export NSPHS_INF=/data/stefane/NSPHS_INF
  ls $NSPHS_INF > sumstats/NSPHS_INF.list
  sort -k2,2 inf1.list > inf1.tmp
  ls $NSPHS_INF | \
  sed 's/NSPHS_inf1_//g;s/.txt.gz//g' | \
  tr '_' '\t' | \
  awk '{
    if(NF==3) $1=$1 "_" $2
    print $1,$NF
  }' | \
  sort -k2,2 | \
  join -j2 inf1.tmp - | \
  parallel -j8 --env NSPHS_INF -C' ' '
  gunzip -c $NSPHS_INF/NSPHS_inf1_{3}_{1}.txt.gz | \
  awk -f tryggve/NSPHS.awk | \
  gzip -f > work/NSPHS.{2}.gz
  '
}

function PIVUS_ULSAM() {
  # PIVIUS and ULSAM, 82 and 85 lines, respectively
  cut -f1 inf1_gene > inf1.tmp

  ls /data/stefang/pivus_ulsam/pivus* | \
  xargs -l -x basename | \
  sed 's/pivus.all.//g;s/.20161128.txt.gz//g' | \
  sort | \
  join - inf1.tmp | \
  awk '{print "pivus.all." $1 ".20161128.txt.gz"}' > sumstats/PIVUS.list

  ls /data/stefang/pivus_ulsam/ulsam* | \
  xargs -l -x basename | \
  sed 's/ulsam.all.//g;s/.20161128.txt.gz//g' | \
  sort | \
  join - inf1.tmp | \
  awk '{print "ulsam.all." $1 ".20161128.txt.gz"}' > sumstats/ULSAM.list
}

function ORCADES_VIS() {
  # ORCADES, VIS, both 91 lines
  ls /data/erimac/ORCADES/ | \
  grep INF1 > sumstats/ORCADES.list

  ls /data/erimac/VIS | \
  grep INF1 > sumstats/VIS.list
}

function STABILITY() {
  # STABILITY
  sort -k2,2 inf1.list > inf1.tmp
  export STABILITY=/data/niceri/Stability_INF1
  ls $STABILITY | \
  sed 's/STABILITY_//g;s/.txt.gz//g' | \
  awk '
  {
    split($1,a,"_")
    if(a[2]=="IL-1" && a[3]="alpha") a[2]=a[2] "_" a[3]
    if(a[2]=="LAP" && a[3]="TGF-beta-1") a[2]=a[2] "_" a[3]
    if(a[2]=="MIP-1" && a[3]="alpha") a[2]=a[2] "_" a[3]
    if(a[2]=="VEGF" && a[3]="A") a[2]=a[2] "_" a[3]
    print a[1],a[2]
  }' | \
  sort -k1,1 | \
  join -11 -22 - inf1.tmp | \
  parallel -j8 --env STABILITY -C' ' '
    cat \
    $STABILITY/STABILITY_{1}_{2}_chr1.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr2.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr3.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr4.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr5.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr6.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr7.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr8.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr9.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr10.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr11.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr12.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr13.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr14.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr15.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr16.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr17.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr18.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr19.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr20.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr21.txt.gz \
    $STABILITY/STABILITY_{1}_{2}_chr22.txt.gz > work/STABILITY.{3}.gz'

  ls work/STABILITY.* > sumstats/STABILITY.list
}

function STANLEY() {
  # STANLEY, only genotype dosages by chromosomal regions, 10699721 lines
  # ls /data/andmala/STANLEY
  # Now it has INF1 results from PLINK for both lah1 and swe6

  sort -k2,2 inf1.list > inf1.tmp
  awk -vFS="\t" -vOFS="\t" '{print $5,$2,$1}' doc/STANLEY_INF_I_annotation.tsv | \
  sort -k1,1 | \
  join -11 -22 -t$'\t' - inf1.tmp | \
  sort -k2,2n | \
  awk -vFS="\t" -vOFS="\t" '{print $2,$4,$3,$1}' > STANLEY.tsv

  # lah1, No. 3 phenotype has no results

  export STANLEY_lah1=/data/andmala/STANLEY_20180911
  awk -vFS="\t" '{print $1, $2}' STANLEY.tsv | \
  parallel -j10 --env STANLEY_lah1 -C' ' ' 
    cat \
    $STANLEY_lah1/STANLEY_lah1_inf_chr1_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr2_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr3_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr4_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr5_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr6_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr7_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr8_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr9_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr10_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr11_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr12_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr13_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr14_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr15_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr16_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr17_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr18_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr19_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr20_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr21_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_lah1/STANLEY_lah1_inf_chr22_pheno{1}.txt.assoc.dosage.gz > work/STANLEY_lah1-{2}.gz'

  rm -f work/STANLEY_lah1-BDNF.gz
  ls work/STANLEY_lah1* > sumstats/STANLEY.list

  # swe6, No. 3 phenotype has no results
  export STANLEY_swe6=/data/andmala/STANLEY_20180911//swe6_inf
  awk -vFS="\t" '{print $1, $2}' STANLEY.tsv | parallel -j10 --env STANLEY_lah1 -C' ' '
    cat \
    $STANLEY_swe6/STANLEY_swe6_inf_chr1_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr2_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr3_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr4_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr5_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr6_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr7_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr8_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr9_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr10_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr11_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr12_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr13_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr14_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr15_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr16_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr17_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr18_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr19_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr20_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr21_pheno{1}.txt.assoc.dosage.gz \
    $STANLEY_swe6/STANLEY_swe6_inf_chr22_pheno{1}.txt.assoc.dosage.gz > work/STANLEY_swe6-{2}.gz'

  rm -f work/STANLEY_swe6-BDNF.gz
  ls work/STANLEY_swe6* >> sumstats/STANLEY.list
}

function list_all () {
  ## list all direcories
  ls -R \
     /data/adumitr/\
     /data/andbre/\
     /data/andmala/\ 
     /data/anekal/\
     /data/anscho/\
     /data/auditlogs/\
     /data/danzie/\
     /data/darzhe/\
     /data/erimac/\
     /data/gussmi/\
     /data/jampet/\
     /data/lasfol/\
     /data/leopad/\
     /data/niceri/\
     /data/sarber/\
     /data/stefane/\
     /data/stefang/\
     /data/sysadmin/\
     /data/thibou/\
     /data/yangwu/ | grep -v cvd | grep -v CVD > inf1.lsR
}

function md()
{
# ready for format.sh
  if [ ! -d sumstats ]; then mkdir sumstats; fi
  for l in $(ls sumstats/*list); do mkdir sumstats/$(echo $(basename $l)|sed 's/.list//g'); done
}

module load parallel/20190122

$1

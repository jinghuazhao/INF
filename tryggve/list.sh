# 20-3-2019 JHZ

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

function EGCUT() {
  # EGCUT_INF by autosomal, female, male
  ls /data/anekal/EGCUT_INF/ | \
  grep inf > sumstats/EGCUT.list
}

function INTERVAL() {
  # INTERVAL, 92 lines
  ls /data/jampet/upload-20170920/ | \
  grep inf | \
  sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g;s/___/ /g' > sumstats/INTERVAL.list
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
     if ($3=="NTF3") $3="NT3"
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
  sed 's/ /\n/g' KORA/KORA.varlist | \
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
  join -11 -22 - prot.list | \
  grep -v P29459 > sumstats/MadCam.list
}

function NSPHS() {
  # NSPHS
  export NSPHS=/data/stefane/NSPHS_INF
  sort -k2,2 inf1.list > inf1.tmp
  ls $NSPHS | \
  sed 's/NSPHS_inf1_//g;s/.txt.gz//g' | \
  tr '_' '\t' | \
  awk '{
    if(NF==3) $1=$1 "_" $2
    print $1,$NF
  }' | \
  sort -k2,2 | \
  join -j2 inf1.tmp - > sumstats/NSPHS.list
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
  grep inf1 doc/olink.prot.list.txt | \
  sed 's/inf1_//g;s/___/\t/g' | \
  awk -vOFS="\t" '{
    l=tolower($1)
    gsub(/mip.1.alpha/,"ccl3",l)
    gsub(/il.10/,"il10",l)
    gsub(/il10ra/,"il.10ra",l)
    gsub(/il10rb/,"il.10rb",l)
    gsub(/il.2/,"il2",l)
    gsub(/il20/,"il.20",l)
    gsub(/il2rb/,"il.2rb",l)
    gsub(/il22.ra1/,"il.22.ra1",l)
    gsub(/il24/,"il.24",l)
    gsub(/il.33/,"il33",l)
    gsub(/il.4/,"il4",l)
    gsub(/il.5/,"il5",l)
    gsub(/il.6/,"il6",l)
    gsub(/il.7/,"il7",l)
    gsub(/il.8/,"il8",l)
    gsub(/il.13/,"il13",l)
    gsub(/il.18/,"il18",l)
    gsub(/il18r1/,"il.18r1",l)
    gsub(/vegf.a/,"vegfa",l)
    print $1,$2,l
  }' | \
  sort -k3,3> inf1.tmp

  ls /data/erimac/ORCADES/ | \
  grep INF1 | \
  sed 's/ORCADES.INF1.//g;s/_rank.tsv.gz//g' | \
  sort | \
  join -11 -23 - inf1.tmp > sumstats/ORCADES.list

  ls /data/erimac/VIS/ | \
  grep INF1 | \
  sed 's/VIS.INF1.//g;s/_rank.tsv.gz//g' | \
  sort | \
  join -11 -23 - inf1.tmp > sumstats/VIS.list
}

function RECOMBINE() {
  export wd=/data/jinhua/data/RECOMBINE
  # tar xfz /data/asahed/RECOMBINE_INF1_pQTLs.tar.gz
  # ls $wd/RECOMBINE_pQTLs__meta_scallop/*gz | \
  # tar xfz /data/asahed/RECOMBINE_INF1_pQTLs_updated_13thMarch_19.tar.gz
  ls $wd/RECOMBINE_INF1_pQTLs_updated_13thMarch_19/*gz | \
  xargs -x -l basename | \
  sed 's/_RECOMBINE.txt.gz//g;s/___/ /g;s/_/ /g' | \
  cut -d' ' -f1-3 | \
  sort -k1,1 -k2,2 -k3,3 | \
  uniq > sumstats/RECOMBINE.list
}

function STABILITY() {
  # STABILITY
  cut -d',' -f2-4 doc/stability_inf_n_nmissing.csv | \
  sed '1d;s/inf_//g;s|\"||g;s/,/ /g;
     s/CCL3/MIP.1.alpha/;
     s/IL10/IL.10/;
     s/IL13/IL.13/;
     s/IL18/IL.18/;
     s/IL33/IL.33/;
     s/IL4/IL.4/;
     s/IL5/IL.5/;
     s/IL6/IL.6/;
     s/IL7/IL.7/;
     s/IL8/IL.8/;
     s/VEGFA/VEGF.A/' | \
  sort -k1,1 | \
  awk '{print $1,$2}' > STABILITY.N
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
  uniq > sumstats/STABILITY.list
}

function STANLEY() {
  # STANLEY, only genotype dosages by chromosomal regions, 10699721 lines
  # ls /data/andmala/STANLEY
  # Now it has INF1 results from PLINK for both lah1 and swe6
  # for both lah1 and swe6, No. 3 phenotype (BDNF) has no results

  sort -k2,2 inf1.list > inf1.tmp
  awk -vFS="\t" -vOFS="\t" '{print $5,$2,$1}' doc/STANLEY_INF_I_annotation.tsv | \
  sort -k1,1 | \
  join -11 -22 -t$'\t' - inf1.tmp | \
  sort -k2,2n | \
  awk -vFS="\t" '!/BDNF/{print $2,$4,$3,$1}' > sumstats/STANLEY.list
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

$1

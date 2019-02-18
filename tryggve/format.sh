# 18-2-2019 JHZ

module load parallel/20170822
export threads=8

# EGCUT_INF -- SNPID has prefix esv for non-rsids
sort -k2,2 inf1.list > inf1.tmp
cat sumstats/EGCUT_INF.list | \
sed 's/EGCUT_autosomal_/_autosomal\t/g;s/EGCUT_X_female_/_X_female\t/g;s/EGCUT_X_male_/_X_male\t/g;s/_inf_280918.txt.gz//g' | \
sort -k2,2 | \
join -j2 - inf1.tmp | \
sort -k2,2 | \
parallel -j$threads -C' ' '
  gunzip -c /data/anekal/EGCUT_INF/EGCUT{2}_{1}_inf_280918.txt.gz | \
  awk "{if(NR>1&&(index(\$1,\"esv\")||index(\$1,\"ss\"))) \$1=\"chr\" \$2 \":\" \$3;print}" | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/EGCUT_INF/EGCUT{2}.{3}.gz'

# INTERVAL
cat $HOME/INF/sumstats/INTERVAL.list | \
sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g;s/___/ /g' | \
parallel -j$threads -C' ' '
   /usr/bin/gunzip -c /data/jampet/upload-20170920/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
   awk -f tryggve/INTERVAL.awk | \
   awk -f tryggve/order.awk | \
   gzip -f > sumstats/INTERVAL/INTERVAL.{1}.gz'

# LifeLinesDeep -- SNPID has no "chr" prefix for non-rsids
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

# KORA

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
cat sumstats/KORA.list | \
parallel -j5 -C' ' '
   zgrep -v Completed KORA/snptest.{1}.out.gz | \
   awk -f tryggve/KORA.awk | \
   gzip -f > sumstats/KORA/KORA.{3}.gz'

# PS -- only 88 proteins without the following,
# BDNF P23560 BDNF
# LIF P15018 LIF
# NRTN Q99748 NRTN
# NTF3 P20783 NT.3

# NSPHS_INF
ls work/NSPHS*gz | \
sed 's|work/NSPHS\.||g;s/\.gz//g' | \
parallel -j$threads --env HOME=$HOME -C' ' '
  gunzip -c work/NSPHS.{}.gz | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/NSPHS_INF/NSPHS.{}.gz'

# PIVUS and ULSAM SNPID has :I/D suffix and VG prefix
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

# ORCADES and VIS
awk -vOFS="\t" '{
  l=tolower($1)
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
}' inf1.list > inf1.tmp

ls /data/erimac/ORCADES/ | \
grep INF1 | \
sed 's/ORCADES.INF1.//g;s/_rank.tsv.gz//g' | \
sort | \
join -a1 -11 -23 - inf1.tmp | \
parallel -j$threads -C' ' '
  gunzip -c /data/erimac/ORCADES/ORCADES.INF1.{1}_rank.tsv.gz | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/ORCADES/ORCADES.{2}.gz'

ls /data/erimac/VIS/ | \
grep INF1 | \
sed 's/VIS.INF1.//g;s/_rank.tsv.gz//g' | \
sort | \
join -a1 -11 -23 - inf1.tmp | \
parallel -j$threads -C' ' '
  gunzip -c /data/erimac/VIS/VIS.INF1.{1}_rank.tsv.gz | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/VIS/VIS.{2}.gz'

# STABILITY
ls work/STABILITY*gz | \
sed 's/work\///g' | \
parallel -j$threads -C' ' '
  gunzip -c work/{} | \
  awk -vOFS="\t" "(NR==1||(\$8>=0.03 && \$8<=0.97 && \$9>=1e-6 && \$11>0 && \$14>=0.4)){print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15}" | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/STABILITY/{}'

# STANLEY_lahl/STANLEY_swe6
ls work/STANLEY*gz | \
sed 's/work\///g' | \
parallel -j$threads -C' ' '
  gunzip -c work/{} | \
  awk -f tryggve/STANLEY.awk | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/STANLEY/{}'

# to pave way for QCGWAS
mkdir $HOME/INF/sumstats/work

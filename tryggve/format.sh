# 13-3-2019 JHZ

export threads=8

# BioFinder
grep -v -w TNF sumstats/BioFinder.list | \
parallel -j5 -C' ' '
   awk -f tryggve/BioFinder.awk /data/andmala/biofinder_inf/rsannot_runGwas_plasmaImp.{1}_zre_INFI.glm.linear | \
   awk -f tryggve/order.awk | \
   gzip -f > sumstats/BioFinder/BioFinder.{3}.gz'
# version with complete data
gunzip -c /data/jinhua/data/BioFinder/rsannot_runGwas_plasmaImp.TNF_zre_INFI.glm.linear.gz | \
awk -f tryggve/BioFinder.awk | \
awk -f tryggve/order.awk | \
gzip -f > sumstats/BioFinder/BioFinder.TNF.gz

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
   sort -k1,1 | \
   join - snpstats/INTERVAL.snpstats | \
   awk -f tryggve/addinfo.awk | \
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
cat sumstats/KORA.list | \
parallel -j5 -C' ' '
   gunzip -c KORA/{1}.gz | \
   awk -f tryggve/KORA.awk | \
   awk -f tryggve/order.awk | \
   gzip -f > sumstats/KORA/KORA.{3}.gz'

# ps -- 89 proteins without the following,
# BDNF P23560 BDNF
# NRTN Q99748 NRTN
# NTF3 P20783 NT.3

# MadCam
cat sumstats/MadCam.list | \
parallel -j5 -C' ' '
   sed 's/CODE_ALLELE_FQ/CODE_ALL_FQ/g' /data/andmala/madcam/MadCAM.{1}.{2}.txt | \
   awk -f tryggve/order.awk | \
   gzip -f > sumstats/MadCam/MadCam.{3}.gz'

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

# RECOMBINE

export rt=/data/jinhua/data/RECOMBINE/RECOMBINE_pQTLs__meta_scallop
sort -k2,2 inf1.list > inf1.tmp
sort -k3,3 sumstats/RECOMBINE.list | \
join -13 -22 - inf1.tmp | \
parallel -j4 --env rt -C' ' '
  (
    for chr in `seq 22`
    do
      gunzip -c $rt/{2}_{3}___{1}_chr${chr}_RECOMBINE.txt.gz | \
      cut -f3-7,9-17
    done
  ) | \
  awk "NR==1||\$1!=SNPID" | \
  awk -vOFS="\t" "
  {
    if (NR>1) 
    {
      CHR=\$2
      POS=\$3
      a1=\$6
      a2=\$7
      if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
      else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
      \$1=snpid
    }
    print
  }" | \
  gzip -f > sumstats/RECOMBINE/RECOMBINE.{4}.gz
'

# STABILITY
ls work/STABILITY*gz | \
sed 's/work\///g' | \
parallel -j$threads -C' ' '
  if [ "{}" == "STABILITY.IFN.gamma.gz" ] || [ "{}" == "STABILITY.IL.22.RA1.gz" ] || [ "{}" == "STABILITY.TSLP.gz" ]; then
     export MAF=0.1
  else
     export MAF=0.03
  fi
  gunzip -c work/{} | \
  awk -vOFS="\t" -vMAF=$MAF -f tryggve/STABILITY.awk | \
  awk -f tryggve/order.awk | \
  gzip -f > sumstats/STABILITY/{}'

for s in lah1 swe6
do
  if [ $s == "lah1" ]; then
     export N=344
  else
     export N=300
  fi
  ls work/STANLEY_${s}*gz | \
  sed 's/work\///g' | \
  parallel -j$threads --env N -C' ' '
    gunzip -c work/{} | \
    awk -vN=$N -f tryggve/STANLEY.awk | \
    awk -f tryggve/order.awk | \
    gzip -f > sumstats/STANLEY/{}'
done

# to pave way for QCGWAS
mkdir $HOME/INF/sumstats/work

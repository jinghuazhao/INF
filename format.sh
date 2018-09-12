# 31-8-2018 JHZ

module load parallel/20170822

if [ ! -d METAL ]; then mkdir METAL; fi

rm -f METAL/METAL.tmp
touch METAL/METAL.tmp
for study in EGCUT INTERVAL LifeLinesDeep NSPHS ORCADES PIVUS ULSAM VIS
do
   ls sumstats/$study | sed 's/'"$study".'//g' | awk -vstudy=$study '{
      s=$1
      gsub(/.gz|@/,"",s);print s " " ENVIRON["HOME"] "/INF/sumstats/" study "/" study "." s ".gz"
   }' >> METAL/METAL.tmp
done
sort -k1,1 METAL/METAL.tmp > METAL/METAL.list
for p in $(cut -f1 inf1.list)
do
   export run=METAL/$p.run
   echo MARKERLABEL SNPID > $run
   echo ALLELELABELS EFFECT_ALLELE REFERENCE_ALLELE >> $run
   echo EFFECTLABEL BETA >> $run
   echo PVALUELABEL PVAL >> $run
   echo WEIGHTLABEL N >> $run
   echo FREQLABEL CODE_ALL_FQ >> $run
   echo STDERRLABEL SE >> $run
   echo SCHEME SAMPLESIZE >> $run
   echo GENOMICCONTROL OFF >> $run
   echo OUTFILE $HOME/INF/METAL/$p- .tbl >> $run
   echo $p | join METAL/METAL.list - | awk '{$1="PROCESS"; print}' >> $run;
   echo ANALYZE >> $run
   echo CLEAR >> $run
done

# NSPHS
sed 's/NSPHS_//g;s/_descriptives.20170411_inclPval.txt.gz//g' sumstats/NSPHS.list | \
parallel -j1 -C' ' 'ln -fs $HOME/INF/work/NSPHS_{}_descriptives.20170411_inclPval.txt.gz sumstats/NSPHS/NSPHS.{}.gz'

# EGCUT and Estonian_biobank
sort -k2,2 inf1.list > inf1.tmp
ls /data/anekal/EGCUT/|grep inf | sed 's/EGCUT_autosomal_//g;s/_inf_140717.txt.gz//g' | grep -v _X | sort | join -11 -22 - inf1.tmp | \
parallel -j1 -C' ' '/usr/bin/ln -sf /data/anekal/EGCUT/EGCUT_autosomal_{1}_inf_140717.txt.gz sumstats/EGCUT/EGCUT.{2}.gz'
ls /data/anekal/Estonian_biobank/|grep inf | sed 's/EGCUT_autosomal_//g;s/_inf_181117.txt.gz//g' | grep -v _X | sort |  join -11 -22 - inf1.tmp | \
parallel -j1 -C' ' '/usr/bin/ln -sf /data/anekal/Estonian_biobank/EGCUT_autosomal_{1}_inf_181117.txt.gz \
   sumstats/Estonian_biobank/Estonian_biobank.{2}.gz'

# PIVUS and ULSAM
ls /data/stefang/pivus_ulsam/pivus* | xargs -l -x basename | sed 's/pivus.all.//g;s/.20161128.txt.gz//g' | sort | \
   join - inf1_gene | cut -d ' ' -f1,3 |\
parallel -j1 -C' ' '/usr/bin/ln -sf /data/stefang/pivus_ulsam/pivus.all.{1}.20161128.txt.gz sumstats/PIVUS/PIVUS.{2}.gz'
ls /data/stefang/pivus_ulsam/ulsam* | xargs -l -x basename | sed 's/ulsam.all.//g;s/.20161128.txt.gz//g' | sort | \
   join - inf1_gene | cut -d ' ' -f1,3 |\
parallel -j1 -C' ' '/usr/bin/ln -sf /data/stefang/pivus_ulsam/ulsam.all.{1}.20161128.txt.gz sumstats/ULSAM/ULSAM.{2}.gz'

# INTERVAL
ls /data/jampet/upload-20170920/ | grep inf | sed 's/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g;s/___/ /g' | \
parallel -j3 -C' ' '/usr/bin/gunzip -c /data/jampet/upload-20170920/INTERVAL_inf1_{1}___{2}_chr_merged.gz | awk -f INTERVAL.awk | gzip -f >\
   sumstats/INTERVAL/INTERVAL.{1}.gz'

# LifeLines
gunzip -c /data/darzhe/LifeLinesDeep.cistranspQTLs.20171220.txt.gz | cut -f1 | awk '(NR>1)' | sed 's/_/ /g' | cut -d' ' -f2 | sort | uniq | \
   join -11 -23  - inf1_gene | parallel -j1 -C' ' 'gunzip -c /data/darzhe/LifeLinesDeep.cistranspQTLs.20171220.txt.gz | \
   awk -vprotein={1} -vFS="\t" -vOFS="\t" "(NR==1||index(\$1,protein))" | gzip -f > sumstats/LifeLinesDeep/LifeLinesDeep.{1}.gz'

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
ls /data/erimac/ORCADES/|grep INF1|sed 's/ORCADES.INF1.//g;s/_rank.tsv.gz//g'| sort | join -a1 -11 -23 - inf1.tmp | \
parallel -j1 -C' ' '/usr/bin/ln -sf /data/erimac/ORCADES/ORCADES.INF1.{1}_rank.tsv.gz sumstats/ORCADES/ORCADES.{2}.gz'
ls /data/erimac/VIS/|grep INF1|sed 's/VIS.INF1.//g;s/_rank.tsv.gz//g'| sort | join -a1 -11 -23 - inf1.tmp | \
parallel -j1 -C' ' '/usr/bin/ln -sf /data/erimac/VIS/VIS.INF1.{1}_rank.tsv.gz sumstats/VIS/VIS.{2}.gz'

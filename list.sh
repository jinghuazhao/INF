# 2-10-2018 JHZ

module load parallel/20170822

if [ ! -d sumstats ]; then mkdir sumstats; fi

## INF list of proteins
grep inf1 doc/olink.prot.list.txt| sed 's/inf1_//g;s/___/\t/g' > inf1.tmp
echo -e "CD6\tP30202\nCD6\tQ8WWJ7\nFGF.5\tP12034" >> inf1.tmp
sort -k1,1 inf1.tmp > inf1.list

## file list
awk -vOFS="\t" '{l=$1;gsub(/\./,"_",$1);print $1,$2,l}' inf1.list > inf1_gene
cut -f1 inf1_gene > inf1.tmp
# NSPHS, tar.gz
mkdir work
cd work
tar xfz /data/stefane/nsphs_20170411.tar.gz
ls *txt | parallel -j4 -C' ' 'gzip -f {}'
cd -
ls work | sed 's/NSPHS_//g;s/_descriptives.20170411_inclPval.txt.gz//g' | sort | join  - inf1.list | \
awk '{print "NSPHS_" $1 "_descriptives.20170411_inclPval.txt.gz"}' > sumstats/NSPHS.list
# PIVIUS and ULSAM, 82 and 85 lines, respectively
ls /data/stefang/pivus_ulsam/pivus* | xargs -l -x basename | sed 's/pivus.all.//g;s/.20161128.txt.gz//g' | \
   sort | join - inf1.tmp | awk '{print "pivus.all." $1 ".20161128.txt.gz"}' > sumstats/PIVUS.list
ls /data/stefang/pivus_ulsam/ulsam* | xargs -l -x basename | sed 's/ulsam.all.//g;s/.20161128.txt.gz//g' | \
   sort | join - inf1.tmp | awk '{print "ulsam.all." $1 ".20161128.txt.gz"}' > sumstats/ULSAM.list
# STABILITY, only CVD2 !
# ls /data/niceri/
# STANLEY, only genotype dosages by chromosomal regions, 10699721 lines
# ls /data/andmala/STANLEY !
# EGCUT, Estonian_biobank, both 18 lines
ls /data/anekal/EGCUT/ | grep inf > sumstats/EGCUT.list
ls /data/anekal/Estonian_biobank|grep inf > sumstats/Estonian_biobank.list
# INTERVAL, 92 lines
ls /data/jampet/upload-20170920/ | grep inf > sumstats/INTERVAL.list
# LifeLines
ls /data/darzhe/LifeLinesDeep.cistranspQTLs.20171220.txt.gz > sumstats/LifeLinesDeep.list
# ORCADES, VIS, both 91 lines
ls /data/erimac/ORCADES/ | grep INF1 > sumstats/ORCADES.list
ls /data/erimac/VIS | grep INF1 > sumstats/VIS.list
## study directory
for l in $(ls sumstats/*list); do mkdir sumstats/$(echo $(basename $l)|sed 's/.list//g'); done

## list all direcories
ls -R\
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

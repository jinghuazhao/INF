# 22-3-2019 JHZ

function METAL_list()
{
# build the lists
  if [ ! -d METAL ]; then mkdir METAL; fi
  (
  for dir in INTERVAL BioFinder EGCUT MadCam KORA NSPHS ORCADES RECOMBINE STABILITY STANLEY VIS
  do
     ls sumstats/$dir | \
     awk -vdir=$dir '{
        s=$1
        gsub(dir,"",s)
        gsub(/EGCUT|NSPHS/,"",s)
        gsub(/\_autosomal|\_X_female|\_X_male|\_lah1.|\_swe6.|.gz|@/,"",s)
        gsub(/^\./,"",s)
        gsub(/@/,"",$1)
        print s " /data/" ENVIRON["USER"] "/INF/sumstats/" dir "/" $1
     }'
  done
  ) | \
  awk '{print $0, NR}' | \
  sort -k1,1 > METAL/METAL.list
}

function METAL_file()
{
# generate METAL command files
  for p in $(cut -f1 inf1.list)
  do
  (
     echo SEPARATOR TAB
     echo COLUMNCOUNTING STRICT
     echo CHROMOSOMELABEL CHR
     echo POSITIONLABEL POS
     echo CUSTOMVARIABLE N
     echo LABEL N as N
     echo TRACKPOSITIONS ON
     echo AVERAGEFREQ ON
     echo MINMAXFREQ ON
     echo ADDFILTER N ">=" 10
     echo ADDFILTER PVAL ">" 0
     echo MARKERLABEL SNPID
     echo ALLELELABELS EFFECT_ALLELE REFERENCE_ALLELE
     echo EFFECTLABEL BETA
     echo PVALUELABEL PVAL
     echo WEIGHTLABEL N
     echo FREQLABEL CODE_ALL_FQ
     echo STDERRLABEL SE
     echo SCHEME STDERR
     echo GENOMICCONTROL OFF
     echo OUTFILE $HOME/INF/METAL/$p- .tbl
     echo $p | \
     join METAL/METAL.list - | \
     sort -k3,3n | \
     awk '{print "PROCESS", $2}'
     echo ANALYZE
     echo CLEAR
  ) > METAL/$p.metal
  done
}

function METAL_QCGWAS()
{
# all in one directory to get ready for QCGWAS
  cd $HOME/INF/sumstats/work
  cat $HOME/INF/METAL/METAL.list | \
  parallel -j10 -C' ' '
    echo {2}; \
    ln -sf {2}
  '
  cd -
}

function METAL_analysis()
{
# conduct the analysis 
# module load metal/20110325 parallel/20170822
  export rt=$HOME/INF
  ls METAL/*.metal | \
  sed 's/.metal//g' | \
  parallel -j5 --env rt -C' ' '
    metal $rt/{}.metal 2>&1 | \
    tee $rt/{}-1.tbl.log; \
    gzip -f $rt/{}-1.tbl
  '
}

function largest_M()
# the union of SNP lists as initially requested by NSPHS
{
  for dir in EGCUT INTERVAL LifeLinesDeep ORCADES PIVUS STABILITY STANLEY ULSAM VIS
  do
    export file=sumstats/$dir/$(ls sumstats/$dir -rS | tail -n 1 | sed 's/@//g;s/.gz//g')
    echo $file
    gunzip -c $file | awk 'NR>1' | cut -f1 | cut -d' ' -f1 > /data/jinhua/M/$dir
  done

  export female=$(ls sumstats/EGCUT_INF/**X_female* -rS | tail -n 1 | sed 's/@//g')
  export male=$(ls sumstats/EGCUT_INF/**X_male* -rS | tail -n 1 | sed 's/@//g')
  gunzip -c $female | awk 'NR>1' | cut -f1 | cut -d' ' -f1 > /data/jinhua/M/EGCUT_X_female
  gunzip -c $male | awk 'NR>1' | cut -f1 | cut -d' ' -f1 > /data/jinhua/M/EGCUT_X_male

  cd /data/jinhua/M
  cat EGCUT_INF EGCUT_X_female EGCUT_X_male INTERVAL ORCADES STABILITY STANLEY VIS | \
  sort | \
  uniq > /data/jinhua/M/M.union
  cd -
}

$1

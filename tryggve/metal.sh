# 4-7-2019 JHZ

function METAL_list()
{
# build the lists
  if [ ! -d METAL ]; then mkdir METAL; fi
  (
  for dir in INTERVAL BioFinder EGCUT KORA NSPHS ORCADES RECOMBINE STABILITY STANLEY VIS
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
  sort -k1,1 | \
  awk '!(/STABILITY/&&/IL.15RA|ST1A1|MCP.3|FGF.5|AXIN1|IL.17A|IL.17C|IL.4|IL.5|IL.10RA|TNF|LIF|IL.13|IL.20RA|IL.24|IL.20|NRTN|ARTN|IL.1.alpha|IL.2RB|IL.33|IFN.gamma|TSLP|IL.22.RA1/)' > METAL/METAL.list
}

function METAL_files()
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
     echo MARKERLABEL SNPID
     echo ALLELELABELS EFFECT_ALLELE REFERENCE_ALLELE
     echo EFFECTLABEL BETA
     echo PVALUELABEL PVAL
     echo WEIGHTLABEL N
     echo FREQLABEL CODE_ALL_FQ
     echo STDERRLABEL SE
     echo SCHEME STDERR
     echo GENOMICCONTROL OFF
     echo LOGPVALUE ON
     echo OUTFILE $HOME/INF/METAL/$p- .tbl
     echo $p | \
     join METAL/METAL.list - | \
     sort -k3,3n | \
     awk '{print "PROCESS", $2}'
     echo ANALYZE HETEROGENEITY
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
  parallel -j4 --env rt -C' ' '
    metal $rt/{}.metal 2>&1 | \
    tee $rt/{}-1.tbl.log; \
    awk "{
       d3=\$13;
       gsub(/?/,\"\",d3)
       if (length(d3 >= 3) && \$18 >= 3500)
          if (\$12 > -9.30103) print;
          else {
             if (\$14 < 30) print;
             else if (d3 == \"nnn\" || d3 == \"ppp\") print
          }
    }" $rt/{}-1.tbl | \
    gzip -f > $rt/{}-1.tbl.gz
    rm $rt/{}-1.tbl
  '
}
# R
# > -log10(5e-10)
# [1] 9.30103
# head -1 METAL/4E.BP1-1.tbl | sed 's|\t|\n|g' | awk '{print "#" NR,$1}'
#1 Chromosome
#2 Position
#3 MarkerName
#4 Allele1
#5 Allele2
#6 Freq1
#7 FreqSE
#8 MinFreq
#9 MaxFreq
#10 Effect
#11 StdErr
#12 log(P)
#13 Direction
#14 HetISq
#15 HetChiSq
#16 HetDf
#17 logHetP
#18 N

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

#!/usr/bin/bash

export dir=~/rds/public_databases/eQTLCatalogue
if [ ! -d ${dir} ];then
   mkdir ${dir}
fi

cd ${dir}
for f in $(sed '1d' ${INF}/eQTLCatalogue/paths.tsv | cut -f4)
do
  echo $f
  echo $f.tbi
  wget ${f}
  wget ${f}.tbi
done
cd -

function region()
(
  grep -f work/INF1.cis work/INF1.merge | \
  awk -vd=1e6 -v OFS='\t' '
  {
    if($3-$2<=2) {$2=$2-d;$3=$3+d}
    if ($2<0) $2=0
    print $1,$2,$3,$5,$6,$8,$9
  }' | \
  sort -k6,6n -k2,2n
)

region | \
bedtools merge | \
wc -l

region | \

curl -X GET http://www.ebi.ac.uk/eqtl/api/chromosomes/4/associations?paginate=False&study=Alasoo_2018&qtl_group=macrophage_naive&quant_method=ge&bp_lower=14737349&bp_upper=16737284

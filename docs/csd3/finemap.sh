# 11-9-2019 JHZ

cd work
tail -n12 *log_sss | \
awk '!/Run/&&/>/{if($1=="==>") {gsub(/.log_sss/," ",$2); printf $2}; $3=$3+0; if($2=="->" && $3>0.1) print $1, $3}' | \
awk '$1=="CCL11-chr3:46250348_C_T"||NF==3' > log_sss.txt
cat log_sss.txt | \
parallel -C' ' 'export pr={1};export k={2}; R --no-save -q <../csd3/finemap.R > {1}-finemap.log'
cd -

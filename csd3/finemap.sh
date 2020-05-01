#!/usr/bin/bash

for j in $(seq 180)
do
  echo ${j}
  sed '1d' work/INF1.merge | cut -f5,6 | awk -v j=${j} 'NR==j'
  csd3/finemap.ini ${j}
done

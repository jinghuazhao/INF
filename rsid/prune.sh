#!/usr/bin/bash

if [ ! -d ${INF}/prune ]; then mkdir ${INF}/prune; fi
for j in $(seq 180)
do
  echo ${j}
  sed '1d' ${INF}/work/INF1.merge-rsid | cut -f5,6 | awk -v j=${j} 'NR==j'
  ${INF}/rsid/prune.inc ${j}
done

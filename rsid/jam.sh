#!/usr/bin/bash

for j in $(seq 180)
do
  echo ${j}
  sed '1d' ${INF}/work/INF1.merge-rsid | cut -f5,6 | awk -v j=${j} 'NR==j'
  ${INF}/rsid/jam.ini ${j}
done

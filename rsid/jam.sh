#!/usr/bin/bash

export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export study=INTERVAL
export data_type=bgen

cd ${INF}/work
cut -d' ' -f1,2 ${sample} | awk 'NR>2' > ${study}.id

module load gcc/6
awk 'NR>1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF --env data_type --env sample -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/jam.ini {3}
'

cd -

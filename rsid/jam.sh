#!/usr/bin/bash

export study=INTERVAL
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export N=4994
export TMPDIR=/tmp

cut -d' ' -f1,2 ${sample} | awk 'NR>2' > ${study}.id
module load gcc/6
awk 'NR>1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF --env sample --env N --env TMPDIR -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/jam.ini {3}
'

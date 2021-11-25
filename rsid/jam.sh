#!/usr/bin/bash

export study=INTERVAL
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export N=4994
export TMPDIR=${HPC_WORK}/work
export dir=${INF}/sentinels

module load gcc/6

cut -d' ' -f1,2 ${sample} | \
awk 'NR>2' > ${INF}/jam/${study}.id
awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF --env study --env sample --env N --env TMPDIR --env dir -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/jam.inc {3}
'
cd ${INF}/jam
echo Missing results
ls *jam.xlsx | awk '{gsub(/-/," ");print $1 "-" $2}' | sort | \
join - <(awk 'NR>1{print $5 "-" $6, NR-1}' ${INF}/work/INF1.merge-rsid | sort) | tr '-' ' '
cd -

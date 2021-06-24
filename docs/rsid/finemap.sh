#!/usr/bin/bash

# At tryggve
# module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.6
# module load bgen/20180807

awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/finemap.inc {3}
'

#!/usr/bin/bash

# At tryggve
# module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.6
# module load bgen/20180807

join -v2 <(sed '1d' sentinels/INF1.jma-rsid | cut -f1,2 | tr '\t' '-' | sort) \
         <(sed '1d' work/INF1.merge-rsid | cut -f5,6 | awk '{print $1 "-" $2,NR}' | sort -k1,1) | \
tr '-' ' ' | \
parallel --env INF -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/finemap.ini {3}
'

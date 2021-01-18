#!/usr/bin/bash

if [ ! -d ${INF}/cs ]; then mkdir ${INF}/cs; fi

awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/cs.inc {3}
'


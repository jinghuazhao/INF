#!/usr/bin/bash

awk 'NR>1 {print $5,$6,NR-1}' work/INF1.merge-rsid | \
parallel --env INF -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/jam.ini {3}
'

#!/bin/bash

source tryggve/analysis.ini

awk 'NR>1' work/INF1_nold.sentinels | \
parallel -j1 -C' ' '
  export p={1}
  export chr={2}
  export pos={3}
  export r={4}
  export start=$(awk -vpos=${pos} "BEGIN{if(pos <= 1e6) print 0; else print pos-1e6}")
  export end=$(awk -vpos=${pos} "BEGIN{print pos+1e6}")
  (
     echo -e "MarkerName\tP-value\tWeight"
     gunzip -c METAL/${p}-1.tbl.gz | \
     awk -vOFS="\t" -vchr=${chr} -vpos=${pos} -vstart=${start} -vend=${end} \
         "(\$1 == chr && \$2 >= start && \$2 < end){split(\$3,a,\"_\");print a[1],10^$12,$\18}" | \
     sort -k1,1 | \
     join -12 -21 work/snp_pos - | \
     awk -vOFS="\t" "{print \$2, \$3, \$4}"
  ) > METAL/${p}-${r}.lz
  cd METAL
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${p}-${r}.lz \
            --plotonly --chr $chr --start $start --end $end --no-date --rundir .
  mv chr${chr}_${start}-${end}.pdf ${p}-${r}.lz.pdf
  pdftopng -r 300 ${p}-${r}.lz.pdf ${p}-${r}
  mv ${p}-${r}-000001.png ${p}-${r}.lz-1.png
  mv ${p}-${r}-000002.png ${p}-${r}.lz-2.png
  cd -
'

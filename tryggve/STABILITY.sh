#!/usr/bin/bash
# STABILITY

module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL

# ls sumstats/STABILITY/*gz | \
# sed 's|sumstats/STABILITY/STABILITY.||g;s/.gz//g' > STABILITY.list
cat STABILITY.list | \
parallel -j4 -C' ' '
  gunzip -c work/STABILITY.{}.gz | \
  awk -vOFS="\t" -vMAF=0 "{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15}" | \
  awk "
  {
    OFS=\"\t\"
    if (NR>1)
    {
      CHR=\$2
      POS=\$3
      a1=\$6
      a2=\$7
      if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
      else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
      \$1=snpid
    }
    print
  }" | \
  gzip -f > sumstats/STABILITY/STABILITY.{}.gz
  export protein={};
  R --no-save -q < STABILITY.R
'
# ::: IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN IL.10RA IL.5 TNF

#!/bin/bash

export job=$1
export TMPDIR=/tmp
export list=${INF}/work/INF1.merge-rsid
export p=$(awk -v job=$job 'NR==job+1{print $5}' ${list})
export r=$(awk -v job=$job 'NR==job+1{print $6}' ${list})
export pr=${p}-${r}
export chr=$(awk -v job=${job} 'NR==job+1{print $8}' ${list})

if [ -f ${INF}/work/${pr}.jma.cojo ]; then
   rm ${INF}/work/${pr}.jma.cojo ${INF}/work/${pr}.ldr.cojo;
fi

function slct_INTERVAL()
{
gcta-1.9 --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr} \
         --cojo-file ${INF}/work/${p}-rsid.ma \
         --chr ${chr} \
         --extract ${pr}.prune \
         --cojo-slct \
         --cojo-p 5e-10 \
         --maf 0.005 \
         --cojo-collinear 0.9 \
         --out ${pr}
}

function slct_INTERVAL_primary()
{
gcta-1.9 --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr} \
         --cojo-file ${INF}/sentinels/${p}-rsid.ma \
         --chr ${chr} \
         --extract ${pr}.prune \
         --cojo-slct \
         --cojo-p 5e-10 \
         --maf 0.005 \
         --cojo-collinear 0.9 \
         --out ${INF}/sentinels/${pr}
}

cd ${INF}/work
slct_INTERVAL_primary
if [ -f ${INF}/sentinels/${pr}.jma.cojo ]; then R --no-save -q < ${INF}/rsid/slct.R > /dev/null; fi
cd -

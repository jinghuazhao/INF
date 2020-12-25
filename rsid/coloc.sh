#!/usr/bin/bash

for r in {1..59}
do
   export r=${r}
   read prot MarkerName < <(awk -vFS="," '$14=="cis"' ${INF}/work/INF1.merge.cis.vs.trans | awk -vFS="," -vr=${r} 'NR==r{print $2,$5}')
   echo ${r} - ${prot} - ${MarkerName}
   export prot=${prot}
   export MarkerName=${MarkerName}
   if [ ! -f ${INF}/coloc/${prot}-${MarkerName}.pdf ] || [ ! -f ${INF}/coloc/${prot}-${MarkerName}.RDS ]; then
     cd ${HPC_WORK}/work
     R --no-save < ${INF}/rsid/coloc.R 2>&1 tee ${prot}-${MarkerName}.log
     ls *tbi | xargs -I {} bash -c "rm {}"
     cd -
   fi
done
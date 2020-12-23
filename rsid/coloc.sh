#!/usr/bin/bash

for r in {1..59}
do
   export r=${r}
   read prot MarkerName < <(awk -vFS="," '$14=="cis"' work/INF1.merge.cis.vs.trans | awk -vFS="," -vr=${r} 'NR==r{print $2,$5}')
   echo ${r} - ${prot} - ${MarkerName}
   export prot=${prot}
   export MarkerName=${MarkerName}
   cd ${HPC_WORK}/work
   R --no-save < ${INF}/rsid/coloc.R 2>&1 tee ${prot}-${MarkerName}.log
   ls *tbi | xargs -I {} bash -c "rm {}"
   cd -
done

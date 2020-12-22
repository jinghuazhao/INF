#!/usr/bin/bash

export M=1e6
for row in {1..180}
do
   read prot MarkerName < <(awk -vrow=${row} 'NR==row+1{print $5,$6}' ${INF}/work/INF1.merge)
   echo ${prot} - ${MarkerName}
   export prot=${prot}
   export MarkerName=${MarkerName}
   cd coloc
   R --no-save < ${INF}/rsid/coloc.R > ${prot}-${MarkerName}.log
   ls *tbi | xargs -I {} bash -c "rm {}"
   cd -
done

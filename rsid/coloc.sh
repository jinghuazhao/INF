#!/usr/bin/bash

for row in {1..59}
do
   export row=${row}
   read prot MarkerName < <(awk -vFS="," '$14=="cis"' work/INF1.merge.cis.vs.trans | awk -vFS="," -vrow=${row} 'NR==row{print $2,$5}')
   echo ${row} - ${prot} - ${MarkerName}
   export prot=${prot}
   export MarkerName=${MarkerName}
   cd ${HPC_WORK}/work
   R --no-save < ${INF}/rsid/coloc.R > ${prot}-${MarkerName}.log
   ls *tbi | xargs -I {} bash -c "rm {}"
   cd -
done

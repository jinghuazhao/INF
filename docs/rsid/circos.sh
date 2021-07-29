#!/usr/bin/bash

export inf_circos=${INF}/circos
export circos_dir=circos-0.69-9
export symdir=${HPC_WORK}/${circos_dir}/inf

if [ ! -d ${inf_circos} ]; then mkdir ${inf_circos}; fi
cd ${inf_circos}

R --no-save -q < ${INF}/rsid/circos.R

if [! -d ${symdir} ]; then mkdir ${symdir}; fi
cd ${symdir}
# for f in $(ls ${inf_circos}/*); do ln -sf $f; done
cp ${inf_circos}/* .
../bin/circos --conf circos.conf


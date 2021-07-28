#!/usr/bin/bash

export inf_circos=${INF}/circos

if [! -d ${inf_circos} ]; then mkdir ${inf_circos}; fi
cd ${inf_circos}

R --no-save -q < ${INF}/rsid/circos.R

export circos_dir=circos-0.69-9
export symdir=${HPC_WORK}/${circos_dir}/circos

if [! -d ${symdir} ]; then mkdir ${symdir}; fi
cd ${symdir}
for f in $(ls $INF/circos/*); do ln -sf $f; done
../bin/circos --conf circos.conf


#!/usr/bin/bash

export circos_dir=circos-0.69-9
export inf_circos=${INF}/circos
export infdir=${HPC_WORK}/${circos_dir}/inf

if [ ! -d ${inf_circos} ]; then mkdir ${inf_circos}; fi
if [ ! -d ${infdir} ]; then mkdir ${infdir}; fi

cd ${inf_circos}
cat ${HPC_WORK}/${circos_dir}/data/karyotype/karyotype.human.hg19.txt | \
awk '/band/&&!/hsX/&&!/hsY/ {gsub(/hs/,"chr",$2);print}' | \
awk -v OFS="\t" '{print $2,$5,$6,$3,$7}' > ${inf_circos}/cytoband.txt

R --no-save -q < ${INF}/rsid/circos2.R

cd ${infdir}
cp ${inf_circos}/* .
sed 's/"//g;s/ //g' ${inf_circos}/pQTL_labels.txt | \
awk -vFS="\t" '{split($4,a,"/");$4=a[1]};1' > pQTL_labels.txt

${circos_dir}/bin/circos --conf circos.conf

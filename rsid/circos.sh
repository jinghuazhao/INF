#!/usr/bin/bash

export inf_circos=${INF}/circos
export circos_dir=circos-0.69-9
export infdir=${HPC_WORK}/${circos_dir}/inf

if [ ! -d ${inf_circos} ]; then mkdir ${inf_circos}; fi
cd ${inf_circos}

awk '/band/&&!/hsX/&&!/hsY/ {gsub(/hs/,"chr",$2);print}' ${HPC_WORK}/${circos_dir}/data/karyotype/karyotype.human.hg19.txt | \
awk -v OFS="\t" '{print $2,$5,$6,$3,$7}' > ${inf_circos}/cytoband.txt

R --no-save -q < ${INF}/rsid/circos.R

convert -density 300 circlize.eps circlize.pdf
convert -density 300 circlize.eps circlize.png

if [ ! -d ${infdir} ]; then mkdir ${infdir}; fi
cd ${infdir}
# for f in $(ls ${inf_circos}/*); do ln -sf $f; done
cp ${inf_circos}/* .
sed 's/"//g' ${inf_circos}/pQTL_labels.txt | \
awk -vFS="\t" '{split($4,a,"/");$4=a[1]};1' > pQTL_labels.txt
../bin/circos --conf circos.conf

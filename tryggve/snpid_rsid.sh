#!/bin/bash

export job=132
export TMPDIR=/data/jinhua/work
export INF=/home/jinhua/INF
export list=${INF}/work/INF1_nold.sentinels
export p=$(awk -v job=${job} 'NR==job+1{print $1}' ${list})
export chr=$(awk -v job=${job} 'NR==job+1{print $2}' ${list})
export pos=$(awk -v job=${job} 'NR==job+1{print $3}' ${list})
export r=$(awk -v job=${job} 'NR==job+1{print $4}' ${list})
export pr=${p}-${r}
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')
export bfile=${INF}/EUR
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export snpid_rsid=${INF}/work/INTERVAL.rsid

cd work
# unpruned
awk -vchr=$chr -vlower=$start -vupper=$end '$1==chr && $4>=lower && $4 < upper{print $2}' ${bfile}.bim | \
sort -k1,1 > ${pr}
join -a2 -e "NA" ${snpid_rsid} ${pr} -o2.1,1.2 > ${pr}.rsid
R --no-save -q <<END
    f <- Sys.getenv("pr")
    snpid_rsid <- read.table(paste0(f,".rsid"), as.is=TRUE, col.names=c("rsid","name"), fill=TRUE)
    nmiss <- with(snpid_rsid,is.na(name))
    snpid_rsid <- within(snpid_rsid, {name[nmiss] <- make.names(rsid[nmiss])})
    save(snpid_rsid, file=paste0(f,".rda"))
END

module load plink/2.00-alpha
# pruned
plink2  --bgen ${bfile}.bgen -sample ${sample} \
        --snp ${r} \
        --window 1000 \
        --maf 0.01 \
        --indep-pairwise 1000kb 1 0.8 \
        --out ${dir}/${pr}
if [ $(grep ${r} ${dir}/${pr}.prune.in | wc -l) -eq 0 ]; then
   (
     echo ${r}
     cat ${dir}/${pr}.prune.in
   ) > $TMPDIR/${pr}
   plink-1.9 --bfile ${bfile} \
             --extract $TMPDIR/${pr} \
             --r2 square \
             --out $TMPDIR/${pr}
   export i=$(paste $TMPDIR/${pr} $TMPDIR/${pr}.ld | cut -f1,2 | sort -r -k2,2g | awk 'NR==1 {print $1}')
   sed -i 's/'"$i"'/'"$r"'/g' ${dir}/${pr}.prune.in
fi
sort -k1,1 ${dir}/${pr}.prune.in > $TMPDIR/${pr}.prune.in
awk 'NR > 1' ${dir}/${p}.ma | \
sort -k1,1 | \
join -j1 - $TMPDIR/${pr}.prune.in | \
awk '{print $1}' > ${dir}/${pr}.prune
join -a2 -e "NA" <(sort -k1,1 ${snpid_rsid}) ${pr}.prune -o2.1,1.2 > ${pr}.prune.rsid
R --no-save -q <<END
    f <- Sys.getenv("pr")
    snpid_rsid <- read.table(paste0(f,".prune.rsid"), as.is=TRUE, col.names=c("rsid","name"), fill=TRUE)
    nmiss <- with(snpid_rsid,is.na(name))
    snpid_rsid <- within(snpid_rsid, {name[nmiss] <- make.names(rsid[nmiss])})
    save(snpid_rsid, file=paste0(f,".prune.rda"))
END
cd -

#!/bin/bash

export job=$1
export list=${INF}/work/INF1.merge-rsid
export p=$(awk -v job=${job} 'NR==job+1{print $5}' ${list})
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
export pr=${p}-${r}
export chr=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $8}')
export pos=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $9}')
export flanking=1e6
export start=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{start=$9-d;if(start<0) start=0; print start}')
export end=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{print $9+d}')
export genotype=${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr}

if [ "${study}" == "" ] || [ "${sample}" == "" ] || [ "${N}" == "" ] || [ "${dir}" == "" ] ; then
   echo Environmental variables such as study, sample, N and dir are not set ... exit
   exit 1
fi

function prune()
{
  module load plink/2.00-alpha
  plink2 --bfile ${genotype} --chr ${chr} --from-bp ${start} --to-bp ${end} \
         --geno 0.1 --mind 0.1 --maf 0.01 --indep-pairwise 1000kb 1 0.1 --out ${pr}-jam
  if [ $(grep -w ${r} ${pr}-jam.prune.in | wc -l) -eq 0 ]; then
     export i=$(grep -w -f ${pr}-jam.prune.in ${genotype}.bim | \
     awk -vpos=${pos} 'function abs(x) {if (x<0) return -x; else return x;}
         {d=$4-pos;print $1, $2, $4, d}' | sort -r -k4,4n | awk 'NR==1 {print $2}')
     sed -i 's/'"$i"'/'"$r"'/g' ${pr}-jam.prune.in
  fi
  if [ ${chr} -eq 19 ]; then
     sort ${pr}-jam.prune.in | join -v1 - NLRP2 > ${pr}-jam.prune
  else
     sort ${pr}-jam.prune.in > ${pr}-jam.prune
  fi
}

function jam()
{
  if [ ! -f ${dir}/${pr}.jma.cojo ]; then
     export k=1
  else
     export k=$(awk "END{print NR-1}" ${dir}/${pr}.jma.cojo)
  fi
  ln -fs ${INF}/finemap/${pr}.z
  awk 'NR>1{print $1}' ${pr}.z | \
  sort | \
  join - ${pr}-jam.prune > ${pr}-jam.incl
  qctool -g ${genotype}.bgen -og ${pr}-jam.bgen -ofiletype bgen -incl-rsids ${pr}-jam.incl
  bgenix -g ${pr}-jam.bgen -index -clobber
  ln -sf ${pr}-jam.bgen.bgi ${pr}-jam.bgi
  R -q --no-save < ${INF}/rsid/jam.R > /dev/null
}

function clean_up()
{
  rm ${pr}-jam.log ${pr}-jam.prune.in ${pr}-jam.prune.out
  rm ${pr}-jam.incl ${pr}-jam.prune
  rm ${pr}-jam.bgen ${pr}-jam.bgen.bgi ${pr}-jam.bgi
}

cd ${INF}/jam
prune
jam
cd -

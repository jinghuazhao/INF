#!/bin/bash

export job=$1
export list=${INF}/work/INF1.merge-rsid
export p=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $5}')
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
export pr=${p}-${r}
export chr=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $8}')
export pos=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $9}')
export flanking=1e6
export start=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{start=$9-d;if(start<0) start=0; print start}')
export end=$(sed '1d' ${list} | awk -v job=${job} -v d=${flanking} 'NR==job{print $9+d}')
export sumstats=$SCALLOP/jp549/olink-merged-output/$(ls $SCALLOP/jp549/olink-merged-output | grep -v cvd | grep ${p}___)
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export study=INTERVAL
export N=4994
export n=1000
export TMPDIR=${HPC_WORK}/work
export dir=${INF}/sentinels
# scaled results
export list=${INF}/work/INF1.merge
export p=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $5}')
export r=$(sed '1d' ${list} | awk -v job=${job} 'NR==job{print $6}')
export sumstats=${INF}/finemapping/${p}-${r}.out
export sample=${INF}/finemapping/${p}-${r}.sample

cd ${INF}/finemap
# z0 - INTERVAL data only excluding NLRP2 region if overlapped and employs P value cutoff to avoid nonsense results
(
  echo MarkerName Chromosome Position Allele2 Allele1 MAF Effect StdErr
# zcat ${sumstats} | \
# awk 'NR > 1' | \
  cat ${sumstats} | \
  awk 'NR > 17' | \
  awk -vchr=${chr} -vstart=${start} -vend=${end} '
  {
    if (($3==chr && $4 >= start && $4 < end && $22 < 1e-3) && !(chr==19 && $4 >= 53296855 && $4 <= 54500000)) {
       chr=$3+0;a1=toupper($5);a2=toupper($6);
       if(a1<a2) snpid=chr ":" $4 "_" a1 "_" a2; else snpid=chr ":" $4 "_" a2 "_" a1;
       if ($2!=".") MarkerName=$2; else MarkerName=snpid;
       if ($24 != "NA" && $25 != "NA") print MarkerName, chr, $4, a2, a1, $19, $24, $25
    }
  } '
) > ${pr}.tmp
awk 'NR > 1{print $1}' ${pr}.tmp > ${pr}.incl

# bgen/bgi - duplicates
qctool -g ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr}.bgen -og ${pr}-${job}.bgen -ofiletype bgen -incl-rsids ${pr}.incl
bgenix -g ${pr}-${job}.bgen -index -clobber
ln -sf ${pr}-${job}.bgen.bgi ${pr}-${job}.bgi
bgenix -g ${pr}-${job}.bgen -list > ${pr}.list
awk 'count[$1]++>0 {print $1}' ${pr}.list | sort | uniq | grep -v -f - ${pr}.tmp > ${pr}.z0
awk 'NR > 1{print $1}' ${pr}.z0 > ${pr}.incl
qctool -g ${pr}-${job}.bgen -og ${pr}.bgen -ofiletype bgen -incl-rsids ${pr}.incl
bgenix -g ${pr}.bgen -index -clobber
ln -sf ${pr}.bgen.bgi ${pr}.bgi
rm ${pr}.tmp ${pr}-${job}.* ${pr}.incl ${pr}.list
# z
(
  join <(awk 'NR > 1' ${pr}.z0 | cut -d' ' -f1,4,5  | sort -k1,1) \
       <(bgenix -g ${pr}.bgen -list 2>&1 | awk 'NR>9 && NF==7'| cut -f2,6,7 | sort -k1,1) | \
  awk '{print $1, ($2!=$4)}'
) > ${pr}.flip
(
  awk 'BEGIN {print "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "flip"}'
  join <(awk 'NR > 1' ${pr}.z0 | sort -k1,1) <(sort -k1,1 ${pr}.flip) | \
  awk '{if($9==1) {t=$4;$4=$5;$5=t};$7=-$7; print}' | awk 'a[$1]++==0'
) > ${pr}.z

if [ ! -f ${dir}/${pr}.jma.cojo ]; then
   export k=1
else
   export k=$(awk "END{print NR-1}" ${dir}/${pr}.jma.cojo)
fi

function ld_finemap()
# ldstore 1.1 and finemap 1.3.1.
{
  (
    echo "z;ld;snp;config;cred;log;n_samples"
    echo "${pr}.z;${pr}.ld;${pr}.snp;${pr}.config;${pr}.cred;${pr}.log;$N"
  ) > ${pr}.master
  ldstore_v1.1 --bcor ${pr}-1 --bgen ${pr}.bgen --n-threads 1
  ldstore_v1.1 --bcor ${pr}-1 --merge 1
  ldstore_v1.1 --bcor ${pr}-1 --matrix ${pr}.ld
  rm ${pr}-1_*
  mv ${pr}-1 ${pr}.bcor
  rm -rf ${pr}.cred* ${pr}.snp ${pr}.config
  finemap_v1.3.1 --sss --in-files ${pr}.master --log --n-causal-snps $k --corr-group 0.7 --group-snps
}

function ld_finemap2()
# ldstore 2 and finemap 1.4
{
  (
    echo "z;ld;bgen;bgi;bcor;bdose;snp;config;cred;log;n_samples"
    echo "${pr}.z;${pr}.ld;${pr}.bgen;${pr}.bgi;${pr}.bcor;${pr}.bdose;${pr}.snp;${pr}.config;${pr}.cred;${pr}.log;$N"
  ) > ${pr}.master2
  ldstore_v2.0_x86_64 --in-files ${pr}.master2 --write-bcor
  ldstore_v2.0_x86_64 --in-files ${pr}.master2 --bcor-to-text
  finemap_v1.4_x86_64 --sss --in-files ${pr}.master2 --log --n-causal-snps $k --corr-group 0.7 --group-snps
}

# xlsx
ld_finemap
gzip -f ${pr}.ld
R -q --no-save < ${INF}/rsid/finemap.R > /dev/null
# rm -f ${pr}.z0 ${pr}.flip ${pr}.master ${pr}.master2 ${pr}.bcor ${pr}.ld.gz ${pr}.bgen ${pr}.bgen.bgi ${pr}.bgi

cd -

#!/bin/bash

#SBATCH --account=PETERS-SL3-CPU
#SBATCH --ntasks=1
#SBATCH --job-name=_fm_INF1
#SBATCH --time=12:00:00
#SBATCH --partition=skylake
#SBATCH --array=1-146%15
#SBATCH --mem=128800
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_fm_%A_%a.out
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_fm_%A_%a.err
#SBATCH --export ALL

export job=${SLURM_ARRAY_TASK_ID}
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export TMPDIR=/tmp
export dir=${INF}/sentinels
export list=${INF}/work/INF1.merge
export p=$(awk -v job=${job} 'NR==job+1{print $5}' ${list})
export chr=$(awk -v job=${job} 'NR==job+1{print $8}' ${list})
export pos=$(awk -v job=${job} 'NR==job+1{print $9}' ${list})
export r=$(awk -v job=${job} 'NR==job+1{print $6}' ${list})
export pr=${p}-${r}
export sumstats=${INF}/METAL/${p}-1.tbl.gz
export flanking=1e6
export start=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{start=pos-flanking;if(start<0) start=0;print start}')
export end=$(awk -vpos=${pos} -vflanking=${flanking} 'BEGIN{print pos+flanking}')
export sample=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export study=INTERVAL
export N=4994

function init()
{
# z0
  (
    zcat ${sumstats} | head -1 | awk '{print $3,$1,$2,$5,$4,"MAF",$10,$11}'
    zcat ${sumstats} | \
    awk 'NR > 1' | \
    sort -k3,3 | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} '
    {
      if ($1==chr && $2 >= start && $2 < end) {
      if ($6 < 0.5) maf = $6; else maf = 1-$6
      if (maf > 0 && maf <= 0.5 && $10 != "NA" && $11 != "NA") print $3, $1, $2, toupper($5), toupper($4), maf, $10, $11
    }
    } ' | \
    join ${pr} -
  ) > ${pr}.z0
  awk 'NR > 1{print $1} ' ${pr}.z0 > ${pr}.incl
# bgen
  qctool -g ${INF}/INTERVAL/INTERVAL-${chr}.bgen -s ${sample} -og ${pr}.bgen -ofiletype bgen -incl-rsids ${pr}.incl
# bgi
  bgenix -g ${pr}.bgen -index -clobber
  ln -sf ${pr}.bgen.bgi ${pr}.bgi
# z
  (
    join <(awk 'NR > 1' ${pr}.z0 | cut -d' ' -f1,4,5 | sort -k1,1) \
         <(bgenix -g ${pr}.bgen -list 2>&1 | awk 'NR>9 && NF==7'| cut -f2,6,7 | sort -k1,1) | \
    awk '{print $1, ($2!=$4)}'
  ) > ${pr}.flip
  (
    awk 'BEGIN {print "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "flip"}'
    join <(awk 'NR > 1' ${pr}.z0) ${pr}.flip | awk '{if($9==1) {t=$4;$4=$5;$5=t};$7=-$7; print}'
  ) > ${pr}.z
}

function gcta()
{
  if [ ! -f ${dir}/${p}.ma ]; then
  (
    echo SNP A1 A2 freq b se p N
    gunzip -c ${sumstats} | \
    awk '{print $3, $4, $5, $6, $10, $11, 10^$12, $18}'
  ) > ${dir}/${p}.ma
  fi
  if [ -f ${dir}/${pr}.jma.cojo ]; then rm ${dir}/${pr}.jma.cojo ${dir}/${pr}.ldr.cojo; fi
  gcta-1.9 --bfile ${INF}/INTERVAL/cardio/INTERVAL \
           --cojo-file ${dir}/${p}.ma \
           --extract ${dir}/${pr}-NLRP2 \
           --cojo-slct \
           --cojo-p 5e-10 \
           --cojo-collinear 0.9 \
           --maf 0.01 \
           --out ${dir}/${pr}
  R -q --no-save < ${INF}/csd3/slct.R > /dev/null
}

function finemap()
{
# master
  (
    echo "z;bgen;bgi;dose;snp;config;cred;log;n_samples"
    echo "${pr}.z;${pr}.bgen;${pr}.bgi;${pr}.dose;${pr}.snp;${pr}.config;${pr}.cred;${pr}.log;$N"
  ) > ${pr}.master
# finemap
  rm -rf ${pr}.cred ${pr}.dose ${pr}.snp ${pr}.config
  if [ ! -f ${dir}/${pr}.jma.cojo ]; then
     export k=1
  else
     export k=$(awk "END{print NR-1}" ${dir}/${pr}.jma.cojo)
  fi
  finemap_v1.3.1 --log --sss --in-files ${pr}.master --n-causal-snps $k --corr-group 0.7 --group-snps
  R -q --no-save < ${INF}/csd3/finemap.R > /dev/null
}

function jam()
{
  export data_type=bgen
  if [ $data_type == "bgen" ]; then cut -d' ' -f1,2 ${sample} | awk 'NR>2' > ${study}.id; fi
  if [ $data_type == "binary_ped" ]; then qctool -g ${pr}.bgen -filetype bgen -s ${sample} -ofiletype binary_ped -og ${pr}; fi
  R -q --no-save < ${INF}/csd3/jam.R > /dev/null
}

cd work
init
gcta
finemap
jam
cd -

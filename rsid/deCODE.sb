#!/usr/bin/bash

#SBATCH --job-name=_deCODE
#SBATCH --account=CARDIO-SL0-CPU
#SBATCH --partition=cardio
#SBATCH --cpus-per-task=6
#SBATCH --qos=cardio
#SBATCH --array=1-72
#SBATCH --mem=10000
#SBATCH --time=3-00:00:00
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_deCODE_%A_%a.o
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_deCODE_%A_%a.e
#SBATCH --export ALL

export TMPDIR=${HPC_WORK}/work
export dir=~/rds/results/public/proteomics/deCODE
export job=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]' ${dir}/urls.txt)
export job=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]' ${INF}/deCODE/olink_inf.lst)

. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load gcc/6
module load aria2-1.33.1-gcc-5.4.0-r36jubs

function downloads()
{
  echo ${job} | aria2c -c -j6 -i -
# aria2c -c -j5 -i pickup.txt
}

function setup()
{
  cd ${dir}/full_downloads
  for f in $(ls ${dir}/sr827_full_download/*md5sum | xargs -l basename -s .md5sum)
  do
    ln -sf ../sr827_full_download/${f}.gz
    ln -sf ../sr827_full_download/${f}.md5sum
  done
  cd -
}

function list()
{
  ls sr827_full_download/*md5sum full_downloads/*md5sum | grep -v "@" | xargs -l basename -s .md5sum | grep -f - -v urls.txt > pickup.txt
  join -13 -22 <(cut -f1,4,7 --output-delimiter=' ' ${INF}/deCODE/SomaLogicv4.tsv | sort -k3,3) \
               <(cut -f2,4,5,20 --output-delimiter=' ' ${INF}/work/INF1.METAL | awk '{print $2":"$3,$4,$1}' | \
  sort -k2,2) | cut -d' ' -f2 | uniq > ${INF}/deCODE/olink_inf.lst
  grep unexpected ${dir}/work/* | tr ':' '\t' | cut -f3 | xargs -l basename | grep -f - ${dir}/urls.txt > ${dir}/191.lst
  grep 3388_58_PAK7_PAK7.txt ${dir}/urls.txt >> ${dir}/191.lst
  ls 190/*md5sum | xargs -l basename -s .md5sum > ${dir}/190.md5sum
  ls 190/*md5sum | xargs -l basename | parallel -C' ' 'diff 190/{} assocs_filtered/{}'

}

function checks()
{
  gunzip -c ${dir}/190/${job}.gz > ${dir}/assocs_filtered/${job}
  md5sum assocs_filtered/${job} > assocs_filtered/${job}.md5sum
  diff ${dir}/190/${job}.md5sum ${dir}/assocs_filtered/${job}.md5sum
  rm ${dir}/assocs_filtered/${job}
}

function olink_inf_idx()
{
# module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7
  cd ${INF}/deCODE
  (
    if [ -f ${dir}/full_downloads/${job}*gz ]; then
       gunzip -c ${dir}/full_downloads/${job}*.gz
    else
       gunzip -c ${dir}/sr827_full_download/${job}*.gz
    fi
  ) | \
  awk "{if(NR>1) gsub(/chr/,\"\",\$1)};1" | \
  tr ' ' '\t' | \
  bgzip -f > ${INF}/deCODE/${job}.gz
  tabix -f -S1 -s1 -b2 -e2 ${INF}/deCODE/${job}.gz
  cd -
}

cd ${dir}
olink_inf_idx
cd -

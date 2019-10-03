# 3-10-2019 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export srcdir=${INF}/ukb

sed '1d' work/INF1.merge | \
sed 's/chr//g' | \
awk -v srcdir=${srcdir} -v flanking=1e6 -v INF=${INF} '
{
  if ($2 >= flanking) start=$2-flanking;
  else start = 0;
  end = $3 + flanking
  range = $1 ":" start "-" end;
  if($1<=9) range0=0 range;
  else range0=range
  cmd=sprintf("qctool -g %s/ukb_imp_chr%d_v3.bgen -incl-range %s -ofiletype bgen -og ukb/%s-chr%s.bgen", srcdir, $1, range0, $5, $6, $5, $6)
  print cmd
}' > work/ukb.list
for i in `seq 22`; do grep chr${i}_ work/ukb.list | wc -l | awk -v chr=${i} '{print chr, $1}'; done
for i in $(seq 22); 
do 
  export i=${i}
  if [ -f ${srcdir}/ukb_imp_chr${i}_v3.bgen ]; then 
     grep chr${i}_ work/ukb.list > work/ukb-${i}.list; 
     export jobs=$(wc -l work/ukb-${i}.list | cut -d' ' -f1)
     (
       echo \#\!/bin/bash
       echo
       echo \#SBATCH --account=PETERS-SL3-CPU
       echo \#SBATCH --ntasks=1
       echo \#SBATCH --job-name=_ukb
       echo \#SBATCH --time=8:00:00
       echo \#SBATCH --cpus-per-task=2
       echo \#SBATCH --partition=skylake
       echo \#SBATCH --mem=128800
       echo \#SBATCH --array=1-${jobs}%10
       echo \#SBATCH --output=/rds/user/jhz22/hpc-work/work/_ukb_%A_%a.out
       echo \#SBATCH --error=/rds/user/jhz22/hpc-work/work/_ukb_%A_%a.err
       echo \#SBATCH --export ALL
       echo
       echo export job=\${SLURM_ARRAY_TASK_ID}
       echo export TMPDIR=/rds/user/jhz22/hpc-work/work
       echo
       echo "awk -v job=\${job} 'NR==job' work/ukb-\${i}.list | bash"
     ) > work/ukb-${i}.sb
     sbatch work/ukb-${i}.sb
  fi; 
done

function combined ()
{
  sed '1d' work/INF1.merge | \
  sortBed -i | \
  mergeBed -i - -d 1000000 | \
  sed 's/chr//g' | \
  awk -v srcdir=${srcdir} -v INF=${INF} '
  {
    range = $1 ":" $2 "-" $3;
    if($1<=9) range0=0 range;
    else range0=range
    cmd=sprintf("qctool -g %s/ukb_imp_chr%d_v3.bgen -incl-range %s -ofiletype bgen -og ukb/chr%s.bgen", srcdir, $1, range0, range)
    print cmd
  }'
}

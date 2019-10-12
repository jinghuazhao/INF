# 12-10-2019 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export intervaldir=${INF}/INTERVAL/per_chr

sed '1d' work/INF1.merge | \
sed 's/chr//g' | \
awk -v intervaldir=${intervaldir} -v flanking=1e6 -v INF=${INF} '
{
  if ($2 >= flanking) start=$2-flanking;
  else start = 0;
  end = $3 + flanking
  range = $1 ":" start "-" end;
  if($1<=9) range0=0 range;
  else range0=range
  cmd=sprintf("qctool -g %s/interval.imputed.olink.chr_%d.bgen -incl-range %s -ofiletype bgen -og INTERVAL/bgen/%s-chr%s.bgen", 
              intervaldir, $1, range0, $5, $6, $5, $6)
  print cmd
}' > work/interval.list
for i in `seq 22`; do grep chr_${i} work/interval.list | wc -l | awk -v chr=${i} '{print chr, $1}'; done
for i in $(seq 22); 
do 
  export i=${i}
  if [ -f ${intervaldir}/interval.imputed.olink.chr_${i}.bgen ]; then 
     grep chr_${i} work/interval.list > work/interval-${i}.list; 
     export jobs=$(wc -l work/interval-${i}.list | cut -d' ' -f1)
     (
       echo \#\!/bin/bash
       echo
       echo \#SBATCH --account=PETERS-SL3-CPU
       echo \#SBATCH --ntasks=1
       echo \#SBATCH --job-name=_INTERVAL
       echo \#SBATCH --time=12:00:00
       echo \#SBATCH --cpus-per-task=2
       echo \#SBATCH --partition=skylake
       echo \#SBATCH --mem=128800
       echo \#SBATCH --array=1-${jobs}%10
       echo \#SBATCH --output=/rds/user/jhz22/hpc-work/work/_INTERVAL_%A_%a.out
       echo \#SBATCH --error=/rds/user/jhz22/hpc-work/work/_INTERVAL_%A_%a.err
       echo \#SBATCH --export ALL
       echo
       echo export job=\${SLURM_ARRAY_TASK_ID}
       echo export TMPDIR=/rds/user/jhz22/hpc-work/work
       echo
       echo "awk -v job=\${job} 'NR==job' work/interval-\${i}.list | bash"
     ) > work/interval-${i}.sb
     sbatch work/interval-${i}.sb
  fi; 
done

sbatch csd3/INTERVAL.sb

function excl()
{
  for chr in 3 6
  do
    bgenix -g ukb/ukb_imp_chr${chr}_v3.bgen -list -incl-range ukb/ukb.range 2>&1 | \
    awk 'NR>9 && NF==7'| \
    cut -f3,4,6,7 | \
    awk '
    {
      CHR=$1
      POS=$2
      a1=$3
      a2=$4
      if (a1>a2) {t=a1; a1=a2; a2=t}
      snpid="chr" CHR+0 ":" POS "_" a1 "_" a2
      print snpid
    }' > ukb/ukb-${chr}.excl
  done
}

function checklist()
{
  cd ukb/nodup
  ls *bgen | sed 's/.bgen//g' > ../../nodup.list
  cd -
  cd ukb/bgen
  ls *bgen | sed 's/.bgen//g' > ../../bgen.list
  cd -
  sdiff bgen.list nodup.list
  sdiff bgen.list nodup.list | awk '/</' | cut -f1 > 112
  awk 'NR>1{print $5 "-" $6}' work/INF1.merge | grep -n -f 112 - | awk '{split($1,a,":");printf a[1] ","}' > 112.list
}

function sentinels_combined ()
{
  sed '1d' work/INF1.merge | \
  sortBed -i | \
  mergeBed -i - -d 1000000 | \
  sed 's/chr//g' | \
  awk -v ukbdir=${ukbdir} '
  {
    range = $1 ":" $2 "-" $3;
    if($1<=9) range0=0 range;
    else range0=range
    cmd=sprintf("qctool -g %s/ukb_imp_chr%d_v3.bgen -incl-range %s -ofiletype bgen -og ukb/chr%s.bgen", ukbdir, $1, range0, range)
    print cmd
  }' > ${INF}/work/ukb.list
}

#!/bin/bash

export TMPDIR=/rds/user/jhz22/hpc-work/work
export rt=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/INTERVAL
export s=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample

cut -d' ' -f1-2,5-7 ${s} | awk 'NR>2' > ${rt}.covar
cut -d' ' -f1-2,4,8-28 ${s} | awk 'NR>2' > ${rt}.qcovar
cut -d' ' -f1-2,29- ${s} | awk 'NR>2' > ${rt}.pheno

for i in $(seq 22);
do
  export i=${i}
  (
    echo \#\!/bin/bash
    echo
    echo \#SBATCH --account=PETERS-SL3-CPU
    echo \#SBATCH --ntasks=1
    echo \#SBATCH --job-name=_ukb
    echo \#SBATCH --time=12:00:00
    echo \#SBATCH --cpus-per-task=2
    echo \#SBATCH --partition=skylake
    echo \#SBATCH --mem=128800
    echo \#SBATCH --array=1-${jobs}%10
    echo \#SBATCH --output=/rds/user/jhz22/hpc-work/work/_h2_%A_%a.out
    echo \#SBATCH --error=/rds/user/jhz22/hpc-work/work/_h2_%A_%a.err
    echo \#SBATCH --export ALL
    echo
    echo export job=\${SLURM_ARRAY_TASK_ID}
    echo export TMPDIR=/rds/user/jhz22/hpc-work/work
    echo
    echo "plink --bgen ${rt}-${i}.bgen \
                --exclude range tryggve/high-LD-regions-hg19.txt \
                --indep-pairwise 500kb 1 0.80 --maf 0.01 --out $rt;\
          plink --bfile ${rt}-${i} --extract ${rt}-${i}.prune.in --make-bed --out ${rt}-${i}.prune" | \
    sed 's/                //g'
  ) > INTERVAL/${rt}-${i}.sb
  sbatch --dry-run INTERVAL/${rt}-${i}.sb
done

seq 22 | awk '{print "INTERVAL/INTERVAL-" $1 ".prune"} > INTERVAL/INTERVAL.list
plink --merge-list INTERVAL/INTERVAL.list --make-bed --out INTERVAL.prune

plink --bfile ${rt}.prune --make-grm-bin --threads 2 --out ${rt}

cut -d' ' -f29- ${s} | head -1 | sed 's/ /\n/g' | awk '{split($1,a,"__"); print a[1]}' > prot.list

sbatch --wait csd3/h2.sb

cd work
grep V\(G\) *hsq | grep Vp | sed 's|.hsq:V(G)/Vp||g' > h2.tsv
cd -

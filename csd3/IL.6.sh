# 24-9-2019 JHZ

export p=IL.6
export srcdir=/rds-d4/user/jhz22/hpc-work/data/ukb
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF

cd ${INF}
#
# exceedingly slow?
# qctool -g ${srcdir}/ukb_imp_chr1_v3.bgen -incl-range 01:153426970-155426970 -ofiletype bgen -og ${p}.bgen
#
module load plink/2.00-alpha
plink2 --bgen ${srcdir}/ukb_imp_chr1_v3.bgen --sample ${srcdir}/ukb20480_imp_chr1_v3_s487395.sample \
       --chr 01 --from-bp 153426970 --to-bp 155426970 \
       --make-bed --out ${p}
plink2 --bfile ${p} --rm-dup force-first list --out dup
plink2 --bfile ${p} --exclude dup.rmdup.list --make-bed --out ${p}_nodup
awk '
{
   CHR=$1
   POS=$4
   a1=$5
   a2=$6
   if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
   else snpid="chr" CHR ":" POS "_" a1 "_" a2
   print snpid, $2
}' ${p}_nodup.bim > ${p}.id

plink --bfile ${p}_nodup --update-name ${p}.id 1 2 --make-bed --out ${p}_snpid

cd -

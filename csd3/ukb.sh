# 24-9-2019 JHZ

export p=IL.6
export srcdir=/rds-d4/user/jhz22/hpc-work/data/ukb
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF

module load plink/2.00-alpha
plink2 --bgen ${srcdir}/ukb_imp_chr1_v3.bgen --sample ${srcdir}/ukb20480_imp_chr1_v3_s487395.sample \
       --chr 01 --from-bp 153426970 --to-bp 155426970 \
       --make-bed --out ${INF}/${p}

awk '
{
   CHR=$1
   POS=$4
   a1=$5
   a2=$6
   if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
   else snpid="chr" CHR ":" POS "_" a1 "_" a2
   print snpid, $2
}' ${INF}/${p}.bim > ${INF}/${p}.id

plink --bfile ${INF}/${p} --update-name ${INF}/${p}.id 1 2 --make-bed --out ${INF}/${p}_snpid

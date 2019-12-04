# 4-12-2019 JHZ

export p=IL.18R1
export r=chr2:102992675_C_T
export pr=${p}-${r}

function gcta()
{
  awk -v pr=${pr} '
  {
    for(i=1;i<=NF;i++)
    ldtable[NR,i]=$i
  }
  END {
    for(i=1;i<=NR;i++)
    {
      if(i==1) for(j=1;j<=NF;j++) printf OFS ldtable[1,j];
      else for(j=1;j<=i;j++) printf OFS ldtable[i,j]
      printf "\n"
    }
  }' work/${pr}.ldr.cojo
}

function ldprint()
{
  awk -v pr=${pr} -v option=${2} '
  {
    if (option=="plink") v=$7;
    if (option=="dosage") v=$11;
    else if (option=="ldstore") v=$15
    ldtable[$3,$6]=ldtable[$6,$3]=v
    ldtable[$3,$3]=ldtable[$6,$6]=1
  }
  END {
    infile=sprintf("%s.snps",pr)
    k=0;while (getline snpid < infile) {k++;snps[k]=snpid}; close(infile)
    printf "SNP"; for(i=1;i<=k;i++) printf OFS snps[i]; printf "\n"
    for(i=1;i<=k;i++)
    {
      printf snps[i]; for(j=1;j<=i;j++) {v=ldtable[snps[i],snps[j]];if (v=="") v="NA"; printf OFS v;} printf "\n"
    }
  }' ${pr}.${1}
}

function ldcalc()
{
  cut -f1 work/${pr}.ldr.cojo | \
  sed '1d' > ${pr}.snps

  (
    head -1 INTERVAL/${p}.ma
    grep -f ${pr}.snps -w INTERVAL/${p}.ma
  ) > ${p}.dat

  (
    echo RSID position chromosome A_allele B_allele
    sed 's/chr//;s/:/ /;s/_/ /g' ${pr}.snps | \
    awk '{print "chr" $1 ":" $2 "_" $3 "_" $4, $2, $1, $3, $4}'
  ) > ${pr}.incl

  ldstore_v1.1 --bcor ${pr}-1 --bgen INTERVAL/nodup/${pr}.bgen --ld-n-samples-avail-prop 0.95 --ld-thold 0.01 --n-threads 1
  ldstore_v1.1 --bcor ${pr}-1 --merge 1
  ldstore_v1.1 --bcor ${pr}-1 --matrix ${pr}.table --incl-variants ${pr}.incl
  rm ${pr}-1_*
  mv ${pr}-1 ${pr}.bcor

R --no-save -q <<END
  pr <- Sys.getenv("pr")
  ldr <- read.table(paste0("work/",pr,".ldr.cojo"),as.is=TRUE,header=TRUE)[,-1]
  ldr
  p <- Sys.getenv("p")
  d <- read.table(paste0(p,".dat"),as.is=TRUE,header=TRUE)
END

  plink --bfile INTERVAL/nodup/${pr} --extract ${pr}.snps --r --out ${pr}

  module load plink/2.00-alpha
  for pair in $(awk '{print $3 "+" $6}' ${pr}.ld)
  do
    plink2 --bgen INTERVAL/nodup/${pr}.bgen --sample INTERVAL/o5000-inf1-outlier_in-r2.sample --extract ${pr}.snps \
           --ld $(echo $pair  | sed 's/+/ / ' ) dosage --out ${pair}
  done

  ldstore_v1.1 --bcor work/INTERVAL/INTERVAL-k-INTERVAL/${pr}.bcor --incl-variants ${pr}.incl --matrix L

  (
    echo title
    for pair in $(awk '{print $3 "+" $6}' ${pr}.ld)
    do
      echo $pair $(grep 'D''' ${pair}.log) | \
      awk '{$4=sqrt($4);print}'
    done
  ) | \
  sed '2d;s/\^2//g' | \
  paste ${pr}.ld - > ${pr}.plink

  (
    echo title
    for pair in $(awk '{print $3 "+" $6}' ${pr}.ld)
    do
       echo $(echo $pair  | sed 's/+/ / ' ) | parallel -C ' ' 'awk -v v1={1} -v v2={2} "\$3 == v1 && \$6 == v2" L'
    done
  ) | \
  paste ${pr}.ld - >  ${pr}.ldstore
}

echo GCTA
gcta
echo PLINK
ldprint plink plink
echo dosage
ldprint plink dosage
echo LDSTORE
ldprint ldstore ldstore

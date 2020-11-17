#/usr/bin/bash

module load gcc/6
function collect()
{ 
  echo ${prefix} -- ${id} -- ${trait}
  (
    cat ${prefix}*result.txt | head -1
    grep -w ${id} ${prefix}*result.txt | grep "Wald ratio"
  ) > ${prefix}-${id}.result
  (
    cat ${prefix}*single.txt | head -1
    grep -w ${id} ${prefix}*single.txt | grep -v -e Egger -e Inverse
  ) > ${prefix}-${id}.single
}
if [ ! -d work/mr/pQTLs ]; then mkdir -p work/mr/pQTLs; fi

cd work/mr/pQTLs
for type in cis trans
do
  export type=${type}
  (
  # rsid prot Allele1 Allele2 Freq1 Effect StdErr log.P. cis.trans
    echo SNP Phenotype effect_allele other_allele eaf beta se pval
    cut -f2,3,6,7,8-11,21 ${INF}/work/INF1.METAL | \
    awk -vtype=${type} 'NR>1 && $9==type {print $1,$2,toupper($3),toupper($4),$5,$6,$7,10^$8}'
  ) > INF1_${type}.ins
  export nrows=$(sed '1d' INF1_${type}.ins | wc -l)
  parallel -C' ' '
  export outcomes={1}
  export row={2}
  R --no-save -q <<\ \ END
    outcomes <- Sys.getenv("outcomes")
    ieugwasr::gwasinfo(id = outcomes)
    row <- Sys.getenv("row")
    type <- Sys.getenv("type")
    ivs <- read.table(paste0("INF1_",type,".ins"),as.is=TRUE,header=TRUE)
    prefix <- paste0("INF1_",outcomes,"-",ivs[row,"Phenotype"],"-",type)
    pQTLtools::pqtlMR(ivs[row,],outcomes,prefix=prefix)
    unlink(paste0(prefix,"-heterogeneity.txt"))
    unlink(paste0(prefix,"-pleiotropy.txt"))
  END
  ' ::: $(cat ${INF}/rsid/mrbase-id.txt) ::: $(seq ${nrows})
  parallel -C' ' '
  export outcomes={1}
  export row={2}
  R --no-save -q <<\ \ END
    outcomes <- Sys.getenv("outcomes")
    row <- Sys.getenv("row")
    type <- Sys.getenv("type")
    ivs <- read.table(paste0("INF1_",type,".ins"),as.is=TRUE,header=TRUE)
    prefix <- paste0("efo_",outcomes,"-",ivs[row,"Phenotype"],"-",type)
    pQTLtools::pqtlMR(ivs[row,],outcomes,prefix=prefix)
    unlink(paste0(prefix,"-heterogeneity.txt"))
    unlink(paste0(prefix,"-pleiotropy.txt"))
  END
  ' ::: $(sed '1d' ${INF}/work/efo.txt | cut -f4) ::: $(seq ${nrows})
  export type=${type}
  export prefix=INF1
  export nrows=$(cat ${INF}/rsid/mrbase-id.txt | wc -l)
  for i in $(seq ${nrows})
  do
    export id=$(awk -vnr=${i} 'NR==nr{print $1}' ${INF}/rsid/mrbase-id.txt)
    collect
  done
  export prefix=efo
  export nrows=$(sed '1d' ${INF}/work/efo.txt | wc -l | cut -d' ' -f1)
  for i in $(seq ${nrows})
  do
    export trait=$(sed '1d' ${INF}/work/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $2}')
    export id=$(sed '1d' ${INF}/work/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $4}')
    collect
  done
done

for prefix in INF1 efo
do
  echo ${prefix}
  export all=$(ls ${prefix}*result.txt | wc -l)
  export p=$(bc -l <<< 0.05/${all})
  awk -vp=${p} -vFS="\t" -vOFS="\t" '$NF<p{split($1,a,"-");print $3,$4,a[5],$6,$7,$8,$9}' ${prefix}*result
done

cd -

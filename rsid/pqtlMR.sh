#/usr/bin/bash

module load gcc/6
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
done
for type in cis trans
do
export type=${type}
export nrows=$(sed '1d' INF1_${type}.ins | wc -l)
parallel -C' ' '
  export outcomes={1}
  export row={2}
  R --no-save <<\ \ END
    INF <- Sys.getenv("INF")
    outcomes <- Sys.getenv("outcomes")
    row <- Sys.getenv("row")
    type <- Sys.getenv("type")
    ivs <- read.table(paste0("INF1_",type,".ins"),as.is=TRUE,header=TRUE)
    prefix <- paste0("INF1_",outcomes,"-",ivs[row,"Phenotype"],"-",type)
    pQTLtools::pqtlMR(ivs[row,],outcomes,prefix=prefix)
  END
' ::: $(cat ${INF}/rsid/mrbase-id.txt) ::: $(seq ${nrows})
parallel -C' ' '
  export outcomes={1}
  export row={2}
  R --no-save <<\ \ END
    INF <- Sys.getenv("INF")
    outcomes <- Sys.getenv("outcomes")
    row <- Sys.getenv("row")
    type <- Sys.getenv("type")
    ivs <- read.table(paste0("INF1_",type,".ins"),as.is=TRUE,header=TRUE)
    prefix <- paste0("efo_",outcomes,"-",ivs[row,"Phenotype"],"-",type)
    pQTLtools::pqtlMR(ivs[row,],outcomes,prefix=prefix)
  END
' ::: $(cut -f4 ${INF}/work/efo.txt) ::: $(seq ${nrows})
done
cd -

#/usr/bin/bash

module load gcc/6
cd work
for type in cis trans
do
  export type=${type}
  (
  # rsid prot Allele1 Allele2 Freq1 Effect StdErr log.P. cis.trans
    echo SNP Phenotype effect_allele other_allele eaf beta se pval
    cut -f2,3,6,7,8-11,21 INF1.METAL | \
    awk -vtype=${type} 'NR>1 && $9==type {print $1,$2,toupper($3),toupper($4),$5,$6,$7,$8}'
  ) > INF1_${type}.ins

  R --no-save <<\ \ END
    INF <- Sys.getenv("INF")
    type <- Sys.getenv("type")
    prefix <- paste0("INF1_",type)
    ivs <- within(read.table(paste0(prefix,".ins"),as.is=TRUE,header=TRUE),{pval=10^pval})
    ids <- scan(paste0(INF,"/rsid/mrbase-id.txt"),what="")
    pQTLtools::pqtlMR(ivs,ids,prefix=prefix)
    efo <- read.delim("efo.txt",as.is=TRUE)
    ids <- with(efo,MRBASEID)
    pQTLtools::pqtlMR(ivs,ids,prefix=paste0("efo_",type))
  END
done
cd -

# cut -f1,2,5,6 --complement INF1_${type}-result.txt | awk -vFS="\t" 'NR==1||$5<0.05'| xsel -i

# 20-3-2019 JHZ

!/BETA/{
  OFS="\t"
  if (NR==2) print  "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  CHR=$1
  POS=$3
  if (substr($2,1,2)=="rs") SNPID=$2; else SNPID="chr" CHR ":" POS
  STRAND="NA"
  EFFECT_ALLELE=$4
  REFERENCE_ALLELE=$5
  CODE_ALL_FQ=$6
  BETA=$8
  SE=$9
  PVAL=$10
  RSQ="NA"
  RSQ_IMP=$7
  IMP="NA"
  if (RSQ_IMP>0.3) print SNPID,CHR,POS,STRAND,N,EFFECT_ALLELE,REFERENCE_ALLELE,CODE_ALL_FQ,BETA,SE,PVAL,RSQ,RSQ_IMP,IMP
}

#1 CHR
#2 SNP, which may have form chrxx_xxx_I/D
#3 BP
#4 A1
#5 A2
#6 FRQ
#7 INFO
#8 BETA
#9 SE
#10 P

# To add STRAND, N and IMP
# map <- read.delim("doc/STANLEY_INF_I_annotation.tsv",as.is=TRUE)
# map[,c(2,5,9,10,13)]
# map <- within(map,id=substr(olink.id,5,))

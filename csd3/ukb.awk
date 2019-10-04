# 4-10-2019 JHZ

{
  alt_ids=$1
  CHR=$3
  SNP=$2
  POS=$4
  a1=$6
  a2=$7
  if (a1>a2) {t=a1; a1=a2; a2=t}
  snpid="chr" CHR+0 ":" POS "_" a1 "_" a2
  print alt_ids, SNP, CHR, POS, $6, $7, SNP, snpid, CHR+0, POS, a1, a2
}

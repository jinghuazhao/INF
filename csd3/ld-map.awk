# 20-8-2019 JHZ

{
  CHR=$1
  SNP=$2
  POS=$4
  a1=$5
  a2=$6
  if (a1>a2) {t=a1; a1=a2; a2=t}
  snpid="chr" CHR ":" POS "_" a1 "_" a2
  print SNP, SNP, CHR, POS, $5, $6, snpid, SNP, CHR, POS, a1, a2
}



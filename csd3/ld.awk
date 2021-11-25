# 2-10-2019 JHZ

{
  if ($1<10) CHR= "0" $1;
  else CHR=$1
  SNP=$2
  POS=$4
  a1=$5
  a2=$6
  if (a1>a2) {t=a1; a1=a2; a2=t}
  snpid="chr" $1 ":" POS "_" a1 "_" a2
  print SNP, SNP, CHR, POS, $5, $6, SNP, snpid, $1, POS, a1, a2
}

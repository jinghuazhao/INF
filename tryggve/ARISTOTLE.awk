# 1-4-2020 JHZ

NR>1 && $9!="NA" && $10!="NA" && $13>0.3 {
  CHR=$2
  POS=$3
  A1=$6
  A2=$7
  if (A1<A2) A1A2= "_" A1 "_" A2; else A1A2="_" A2 "_" A1
  SNPID="chr" CHR ":" POS A1A2
  $1=SNPID
  print
}

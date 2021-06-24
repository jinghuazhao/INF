# 6-11-2019 JHZ

{
  if (NR>1)
  {
    pval=10^$13
    if (pval==0 && $26>1) pval=1
    print $4, toupper($5), toupper($6), $7, $11, $12, pval, $19
  }
}

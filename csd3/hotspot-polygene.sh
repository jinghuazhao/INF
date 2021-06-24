# 18-10-2019 JHZ

# long format
cut -d' ' -f3,4,5,7 work/INF1.merge.cis.vs.trans | awk -v OFS="\t" 'NR>1{print "chr" $1,$2-1, $2, $3, $4}' | sed 's/\"//g' > a1

(
  sort -k1,1n -k2,2n csd3/glist-hg19 | \
  grep -v X | \
  grep -v Y | \
  awk '{$1="chr" $1;print}' | \
  sed 's/ /\t/g'
) > a2

bedtools intersect -a a1 -b a2 -wa -wb -loj > sentinel-glist.long
rm a1 a2

# compact format
R --no-save -q <<END
  options(width=1000)
  d <- read.table("sentinel-glist.long",as.is=TRUE)
  a <- aggregate(d,by=list(with(d,V1),with(d,V2),with(d,V3)),FUN="c")
  sink("sentinel-glist.wide")
  a
  sink()
  for (v in c(1:4,6:9))
  {
#     s <- lapply(a[v],unique,"[[")
#     a[v] <- unlist(lapply(s,paste,collapse=","))
      s <- lapply(a[,3+v],unique,"[[")
      a[,3+v] <- unlist(lapply(s,paste,collapse=","))
  }
  sink("sentinel-glist.short")
  a[paste0("V",1:9)]
  sink()
END

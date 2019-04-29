# 29-4-2019 JHZ

# code to summarise protein-region results
(
  echo -e "prot\tregion\t$(cat *jma.cojo | head -1)"
  cat INF1.ailist | \
  parallel -C' ' 'awk "NR>1 {print p, r, \$0}" OFS="\t" p={1} r={2} {1}-{2}.jma.cojo'
) > INF1.jma

R --no-save -q <<END
  jma <- read.delim("INF1.jma")
  jma <- within(jma, {r2=LD_r^2})
  dim(jma)
  subset(jma,r2>0.3)
  subset(jma,r2>0.1)
END

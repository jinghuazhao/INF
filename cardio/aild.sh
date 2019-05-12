# 12-5-2019 JHZ

# code to summarise protein-region results
(
  echo -e "prot\tregion\t$(cat *jma.cojo | head -1)"
  cat work/INF1.aild | \
  parallel -j1 -C' ' 'awk "NR>1 {print p, r, \$0}" OFS="\t" p={1} r={2} aild/cojo/{1}-{2}.jma.cojo'
) > INF1.txt

R --no-save -q <<END
  jma <- read.delim("INF1.txt")
  jma <- within(jma, {r2=LD_r^2})
  dim(jma)
  subset(jma,r2>0.3)
  subset(jma,r2>0.1)
END

# for tryggve/analysis.sh fp

cut -f1,3-16 INF1.txt > INF1.jma

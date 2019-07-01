# 1-7-2019 JHZ

(
  cat dist/cojo/*cma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  cut -d' ' -f1,4 dist/INF1_nold.sentinels | \
  parallel -j1 -C' ' '
    if [ -f dist/cojo/{1}-{2}.cma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "\$NF <= 5e-10 {print prot, snpid, \$0}" dist/cojo/{1}-{2}.cma.cojo
    fi'
) > INF1.cma

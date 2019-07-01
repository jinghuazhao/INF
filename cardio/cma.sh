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

R --no-save -q <<END

options(echo=FALSE)
require(gap)
require(reshape)
cma <- read.delim("INF1.cma", as.is=TRUE)
cma <- rename(cma,c(Chr="Chrom",bp="End",bC="Effect",bC_se="StdErr",SNP="MarkerName"))
tag <- "_cma"

for (protein in with(cma, unique(prot)))
{
  p <- subset(cma, prot==protein)
  chrs <- with(p,unique(Chrom))
  for(chr in chrs)
  {
    ps <- subset(p,Chrom==chr)
    row.names(ps) <- 1:nrow(ps)
    sentinels(ps, protein, 1)
  }
}

END

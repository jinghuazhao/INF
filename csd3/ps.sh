# 6-12-2019 JHZ

for catalogue in eQTL pQTL mQTL methQTL GWAS
do
  export catalogue=${catalogue}
  R --no-save -q <csd3/ps.R > ps.log
done

# gene
R --no-save -q < csd3/ps.gene.R > ps.gene.log

phenoscanner -s T -c All -x EUR -p 0.0000001 --r2 0.6 -i INF1.merge.snp -o INF1

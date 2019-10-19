# 19-10-2019 JHZ

for catalogue in eQTL pQTL mQTL methQTL GWAS
do
  export catalogue=${catalogue}
  R --no-save -q <csd3/ps.R > ps.log
done

# gene
R --no-save -q < csd3/ps.gene.R > ps.gene.log

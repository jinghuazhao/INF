# 18-10-2019 JHZ

for catalogue in eQTL pQTL mQTL methQTL GWAS
do
  export catalogue=${catalogue}
  R --no-save -q <csd3/ps.R > ps.log
done

# gene
R --no-save -q <<END
  glist <- scan("work/INF1.merge.gene",what="")
  genes <- data.frame()
  results <- data.frame()
  for(s in 1:6)
  {
    gset <- ifelse(s==6, 61:68, (s-1)*10+1:10)
    g <- phenoscanner::phenoscanner(genequery=glist[gset], catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
    genes <- rbind(genes,with(g,genes))
    results <- rbind(results,with(g,results))
  }
  r <- list(genes=genes,results=results)
  lapply(r,dim)
  save(r,file="work/INF1.merge.genes")
END

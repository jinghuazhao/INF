# 18-10-2019 JHZ

# To set it up
# explicit rsid
(
  join <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1 | uniq) work/INTERVAL.rsid
) > work/INF1.merge.rsid

# Novel loci
R --no-save -q <<END
  require(phenoscanner)
  rsid <- with(read.table("work/INF1.merge.rsid",as.is=TRUE,col.names=c("snpid","rsid")), rsid)
  r1 <- phenoscanner(snpquery=rsid[1:100], catalogue="pQTL", proxies = "EUR", pvalue = 1e-07, r2= 0.8, build=37)
  lapply(r1,dim)
  r2 <- phenoscanner(snpquery=rsid[101:162], catalogue="pQTL", proxies = "EUR", pvalue = 1e-07, r2= 0.8, build=37)
  lapply(r2,dim)
  r <- list(snps=rbind(with(r1,snps),with(r2,snps)),results=rbind(with(r1,results),with(r2,results)))
  lapply(r,dim)
  save(r,file="work/INF1.merge.pQTL")
END

# SH2B3 and chr12:111884608_C_T sentinel
R --no-save -q <<END
  gene <- phenoscanner::phenoscanner(genequery="TNFSF10", catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(gene,dim)
  snpid <- "chr12:111884608"
  cat(rsid,"\n")
  GWAS <- phenoscanner::phenoscanner(snpquery=snpid, catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(GWAS)
  with(GWAS, results)[c("rsid","study","trait")]
  pQTL <- phenoscanner::phenoscanner(snpquery=snpid, catalogue="pQTL", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(pQTL,dim)
  with(pQTL, results)[c("rsid","study","trait")]
END

# IL.12B
R --no-save -q <<END
  gene <- phenoscanner::phenoscanner(genequery="IL12B", catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(gene,dim)
  g <- with(gene,genes)
  r <- with(gene,results)
END

# PD.L1
R --no-save -q <<END
  gene <- phenoscanner::phenoscanner(genequery="CD274", catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(gene,dim)
  g <- with(gene,genes)
  r <- with(gene,results)
END

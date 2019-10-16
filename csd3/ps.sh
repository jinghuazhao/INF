# 16-10-2019 JHZ

# To set it up
# explicit rsid
(
  join <(sed '1d' work/INF1.merge | cut -f6 | sort -k1,1 | uniq) work/INTERVAL.rsid
) > work/INF1.merge.rsid

# All of the sentinels
R --no-save -q <<END
  require(phenoscanner)
  rsid <- with(read.table("work/INF1.merge.rsid",as.is=TRUE,col.names=c("snpid","rsid")), rsid)
  r1 <- phenoscanner(snpquery=rsid[1:100], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(r1,dim)
  r2 <- phenoscanner(snpquery=rsid[101:162], proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(r2,dim)
  r <- list(rbind(with(r1,snps),with(r2,snps)),rbind(with(r1,results),with(r2,results)))
  lapply(r,dim)
  save(r,file="work/INF1.merge.ps")
END

# Only SH2B3 through chr12:111884608_C_T
R --no-save -q <<END
  snpid <- "chr12:111884608"
  cat(rsid,"\n")
  GWAS <- phenoscanner::phenoscanner(snpquery=snpid, catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(GWAS)
  with(GWAS, results)[c("rsid","study","trait")]
  pQTL <- phenoscanner::phenoscanner(snpquery=snpid, catalogue="pQTL", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(pQTL,dim)
  with(pQTL, results)[c("rsid","study","trait")]
END

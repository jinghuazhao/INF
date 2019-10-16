# 16-10-2019 JHZ

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

# only SH2B3 through chr12:111884608_C_T
R --no-save -q <<END
  rsid <- with(subset(read.table("work/INF1.merge.rsid",as.is=TRUE,col.names=c("snpid","rsid")),snpid=="chr12:111884608_C_T"), rsid)
  cat(rsid,"\n")
  SH2B3 <- phenoscanner::phenoscanner(snpquery=rsid, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(SH2B3,dim)
END

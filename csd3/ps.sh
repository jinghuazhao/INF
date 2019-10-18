# 18-10-2019 JHZ

# lookup
R --no-save -q <<END
  require(phenoscanner)
  rsid <- with(read.table("work/INF1.merge.rsid",as.is=TRUE,col.names=c("snpid","rsid")), rsid)
  for (catalogue in c("eQTL","pQTL","mQTL","methQTL","GWAS"))
  {
    r1 <- phenoscanner(snpquery=rsid[1:100], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.8, build=37)
    lapply(r1,dim)
    r2 <- phenoscanner(snpquery=rsid[101:162], catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.8, build=37)
    lapply(r2,dim)
    r <- list(snps=rbind(with(r1,snps),with(r2,snps)),results=rbind(with(r1,results),with(r2,results)))
    lapply(r,dim)
    save(r,file=paste0("work/INF1.merge.",catalogue))
  }
  m <- read.delim("work/INF1.merge",as.is=TRUE)[c("prot","MarkerName")]
  names(m) <- c("prot","ref_snpid")
END

R --no-save -q <<END
  options(width=500)
  for (catalogue in c("eQTL","pQTL","mQTL","methQTL","GWAS"))
  {
    load(paste0("work/INF1.merge.",catalogue))
    attach(r)
    results <- within(results,{
      a1 <- ref_a1
      a2 <- ref_a2
      swap <- ref_a1 > ref_a2
      a1[swap] <- ref_a2[swap]
      a2[swap] <- ref_a1[swap]
      ref_snpid <- paste0(ref_hg19_coordinates,":",a1,"_",a2)
    })
    for(d in unique(with(results,dataset)))
    {
      cat(d,"\n")
      sink(paste(catalogue,d,sep="."))
      s <- subset(results[c("ref_rsid","ref_snpid","rsid","r2","trait","dataset","pmid")],dataset==d)
      print(s)
      sink()
    }
    detach(r)
  }
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

# PD.L1
R --no-save -q <<END
  gene <- phenoscanner::phenoscanner(genequery="CD274", catalogue="GWAS", proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
  lapply(gene,dim)
  g <- with(gene,genes)
  r <- with(gene,results)
END

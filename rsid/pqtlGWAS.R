options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{hg19_coordinates=paste0("chr",Chromosome,":",Position)})

INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,prot,hg19_coordinates,MarkerName,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.) %>%
  group_by(Chromosome,Position,MarkerName,hg19_coordinates,rsid) %>%
  summarise(n=n(), 
            prots=paste(prot, collapse=";"),
            A1=paste(toupper(Allele1), collapse=";"),
            A2=paste(toupper(Allele2), collapse=";"),
            EAF=paste(Freq1, collapse=";"),
            Effects=paste(Effect, collapse=";"),
            SEs=paste(StdErr, collapse=";"),
            log10P=paste(log.P., collapse=";")
  ) %>%
  data.frame()

rsid <- INF1_aggr[["rsid"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
save(INF1_aggr,r,file=file.path(INF,"work","INF1.merge.GWAS"))

sr <- with(r,
{
  names(snps) <- paste0("s.",names(snps))
  names(results) <- paste0("r.",names(results))
  merge(snps,results,by.x="s.hg19_coordinates",by.y="r.hg19_coordinates",all.y=TRUE)
})
isr <- merge(INF1_aggr,sr,by.x="hg19_coordinates",by.y="s.snpid")
write.csv(isr,file="isr.csv",row.names=FALSE,quote=FALSE)

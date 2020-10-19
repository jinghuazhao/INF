options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                     rsID=rsid;
                     hg19_coordinates=paste0("chr",Chromosome,":",Position)})

INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,prot,uniprot,hg19_coordinates,MarkerName,rsID,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
  group_by(Chromosome,Position,MarkerName,rsID,hg19_coordinates) %>%
  summarise(nprots=n(),
            UniProts=paste(uniprot,collapse=";"),
            prots=paste(prot,collapse=";"),
            A1=paste(toupper(Allele1),collapse=";"),
            A2=paste(toupper(Allele2),collapse=";"),
            EAF=paste(Freq1,collapse=";"),
            Effects=paste(Effect,collapse=";"),
            SEs=paste(StdErr,collapse=";"),
            log10P=paste(log.P.,collapse=";"),
            cistrans=paste(cis.trans,collapse=";")
  ) %>%
  data.frame()

rsid <- INF1_aggr[["rsID"]]
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37

for (catalogue in c("eQTL","mQTL","pQTL"))
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps,results))
  f <- paste0(file.path(INF,"work","INF1.merge."),catalogue)
  save(INF1_aggr,r,ps,file=f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  write.table(ips,file=paste0(f,catalogue,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}

# rework
INF <- Sys.getenv("INF")
for (catalogue in c("eQTL","mQTL","pQTL"))
{
  f <- paste0(file.path(INF,"work","INF1.merge."),catalogue)
  load(f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}


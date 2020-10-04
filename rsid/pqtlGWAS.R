options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                     rsID=rsid;
                     hg19_coordinates=paste0("chr",Chromosome,":",Position)})

INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,prot,hg19_coordinates,MarkerName,rsID,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
  group_by(Chromosome,Position,MarkerName,rsID,hg19_coordinates) %>%
  summarise(nprots=n(),
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
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8 / nrow(INF1_metal)
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
ps <- with(r,right_join(snps,results))
save(INF1_aggr,r,ps,file=file.path(INF,"work","INF1.merge.GWAS"))

efo_list_immune <- subset(read.csv("work/efo_list_annotated.csv",as.is=TRUE),immune_mediated==1)
isd1 <- subset(merge(INF1_aggr,subset(ps,efo%in%with(efo_list_immune,EFO)),by="hg19_coordinates"),select=-c(Chromosome,Position))
write.table(isd1,file="isd1.tsv",row.names=FALSE,quote=FALSE,sep="\t")

load("work/efo.rda")
efo_0000540 <- gsub(":","_",as.data.frame(isd)[["efo_0000540"]])
isd2 <- subset(merge(INF1_aggr,subset(ps,efo%in%efo_0000540),by="hg19_coordinates"),select=-c(Chromosome,Position))
write.table(isd2,file="isd2.tsv",row.names=FALSE,quote=FALSE,sep="\t")

isd2[c("MarkerName","rsID","prots","trait","efo","study","pmid","ancestry","year","beta","se")]

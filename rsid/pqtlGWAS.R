options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{hg19_coordinates=paste0("chr",Chromosome,":",Position)})

INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,prot,hg19_coordinates,MarkerName,rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
  group_by(Chromosome,Position,MarkerName,hg19_coordinates,rsid) %>%
  summarise(n=n(), 
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

rsid <- INF1_aggr[["rsid"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
save(INF1_aggr,r,file=file.path(INF,"work","INF1.merge.GWAS"))

efo_list_immune <- read.csv("work/efo_list_annotated.csv",as.is=TRUE)
EFO <- efo_list_immune %>% select(EFO) %>% summarise(EFO=paste(gsub(","," ",EFO),collapse=" "))
sr <- with(r,
{
  snps <- subset(snps,select=-c(ref_hg38_coordinates,ref_pos_hg38,pos_hg38,hg38_coordinates))
  names(snps) <- paste0("s.",names(snps))
  results <- subset(results,select=-c(ref_hg38_coordinates,hg38_coordinates))
  names(results) <- paste0("r.",names(results))
  merge(snps,results,by.x="s.hg19_coordinates",by.y="r.hg19_coordinates",all.y=TRUE)
})
isr <- merge(INF1_aggr,subset(sr,r.efo%in%with(efo_list_immune,EFO)),by.x="hg19_coordinates",by.y="s.hg19_coordinates")
write.table(isr,file="isr.txt",row.names=FALSE,quote=FALSE,sep="\t")
MarkerName_notfound <- setdiff(INF1_aggr[["MarkerName"]],names(table(isr$MarkerName)))
hg19_coordinates_notfound <- setdiff(INF1_aggr[["hg19_coordinates"]],isr$s.ref_hg19_coordinates)
r_notfound <- snpqueries(hg19_coordinates_notfound, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r_notfound,dim)
r_notfound_efo <- with(r_notfound,efo)


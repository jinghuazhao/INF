options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                     hg19_coordinates=paste0("chr",Chromosome,":",Position)}) %>%
              rename(INF1_rsid=rsid) %>%
              rename(Total=N)

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

rsid <- INF1_aggr[["INF1_rsid"]]
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

# loading
INF <- Sys.getenv("INF")
for (catalogue in c("eQTL","mQTL","pQTL"))
{
  f <- paste0(file.path(INF,"work","INF1.merge."),catalogue)
  load(f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),
                select=-c(Chromosome, Position, EAF, Effects, SEs, nprots, snpid,
                          hg19_coordinates,hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                          ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                          rsid, chr, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                          dprime, efo, n, n_studies, unit, direction))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}

SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,entGene,ensGene,extGene) %>% rename(trait=TargetFullName)
pQTL <- dplyr::left_join(ips,SL)
write.table(pQTL,file=paste0(f,"-SomaLogic.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

# a few checks
gap::pvalue(-1.102/0.016)
gap::pvalue(-1.114/0.0157)
rs12075 <-c("P51671","P80162","P13500","P80075","P80098","Q99616")
subset(SomaLogic160410,UniProt%in%rs12075)
subset(SomaLogic160410,TargetFullName%in%c("Eotaxin","Corneodesmosin","C-C motif chemokine 14","C-C motif chemokine 26"))
subset(SomaLogic160410,TargetFullName=="SLAM family member 7")
subset(SomaLogic160410,UniProt=="P50591")
subset(SomaLogic160410,UniProt=="O14625")
subset(SomaLogic160410,UniProt=="P15692")

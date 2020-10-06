options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                     hg19_coordinates=paste0("chr",Chromosome,":",Position)}) %>%
              rename(INF1_rsid=rsid) %>%
              rename(Total=N)

INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,prot,hg19_coordinates,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans,INF1_rsid) %>%
  group_by(Chromosome,Position,MarkerName,INF1_rsid,hg19_coordinates) %>%
  summarise(nprots=n(),
            prots=paste(prot,collapse=";"),
            Allele1=paste(toupper(Allele1),collapse=";"),
            Allele2=paste(toupper(Allele2),collapse=";"),
            EAF=paste(Freq1,collapse=";"),
            Effects=paste(Effect,collapse=";"),
            SEs=paste(StdErr,collapse=";"),
            log10P=paste(log.P.,collapse=";"),
            cistrans=paste(cis.trans,collapse=";")
  ) %>%
  data.frame()

rsid <- INF1_aggr[["INF1_rsid"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8 / nrow(INF1_metal)
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
ps <- subset(with(r,right_join(snps,results)),select=-c(hg38_coordinates,ref_hg38_coordinates,pos_hg38,ref_pos_hg38,dprime))
save(INF1_aggr,r,ps,file=file.path(INF,"work","INF1.merge.GWAS"))

isd123 <- function()
{
  efo_list_immune <- subset(read.csv("work/efo_list_annotated.csv",as.is=TRUE),immune_mediated==1)
  isd1 <- merge(aggr,subset(ps,efo%in%with(efo_list_immune,EFO)),by="hg19_coordinates")
  write.table(isd1,file="isd1.tsv",row.names=FALSE,quote=FALSE,sep="\t")

  load("work/efo.rda")
  efo_0000540 <- gsub(":","_",as.data.frame(isd)[["efo_0000540"]])
  isd2 <- merge(aggr,subset(ps,efo%in%efo_0000540),by="hg19_coordinates")
  write.table(isd2,file="isd2.tsv",row.names=FALSE,quote=FALSE,sep="\t")

  fang_efo <- gsub(":","_",with(read.delim("doc/fang.efos.txt",as.is=TRUE),id))
  isd3 <- merge(aggr,subset(ps,efo%in%fang_efo),by="hg19_coordinates")
  write.table(isd3,file="isd3.tsv",row.names=FALSE,quote=FALSE,sep="\t")
}

metal <- subset(within(INF1_metal,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
                select=-c(Chromosome,Position,INF1_rsid,Direction))
aggr <- subset(within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
               select=-c(Chromosome,Position,INF1_rsid))
jimmy <- read.delim("doc/immune.efos.txt",as.is=TRUE)
jimmy_efo <- gsub(":","_",with(jimmy,id))
long <- merge(metal,subset(ps,efo%in%jimmy_efo),by="hg19_coordinates")
short <- merge(aggr,subset(ps,efo%in%jimmy_efo),by="hg19_coordinates")

require(openxlsx)
xlsx <- "work/pqtlgwas.xlsx"
wb <- createWorkbook(xlsx)
addWorksheet(wb, "METAL")
writeDataTable(wb, "METAL", subset(INF1_metal,select=-c(INF1_rsid,hg19_coordinates)))
addWorksheet(wb, "ps")
writeDataTable(wb, "ps", ps)
addWorksheet(wb, "EFO")
writeDataTable(wb, "EFO", jimmy)
addWorksheet(wb, "long")
writeDataTable(wb, "long", 
               subset(long,select=-c(hg19_coordinates,ref_hg19_coordinates,ref_protein_position,ref_amino_acids,snp,snpid,protein_position,amino_acids)))
addWorksheet(wb, "short")
writeDataTable(wb, "short", 
               subset(short,select=-c(hg19_coordinates,ref_hg19_coordinates,ref_protein_position,ref_amino_acids,snp,snpid,protein_position,amino_acids)))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

view <- function(id,efoid,
                 v=c("MarkerName","Allele1","Allele2","a1","a2","efo","ref_a1","ref_a2","proxy","r2","beta","se","p","trait","ancestry","pmid","study"))
{
  options(width=200)
  cat(id,efoid,"\n")
  d <- subset(short[v],MarkerName==id & efo==efoid)
  subset(d,select=-c(MarkerName,efo))
}
# celiac
view("chr12:111884608_C_T", "EFO_0001060")
# T1D
view("chr12:111884608_C_T", "EFO_0001359")
# Hypothyroidism
view("chr12:111884608_C_T", "EFO_0004705")
# Allergic disease asthma hay fever or eczema
view("chr12:111932800_C_T", "EFO_0003785")
# Celiac
view("chr12:112007756_C_T","EFO_0001060")
# Crohn's
view("chr19:49206172_C_T","EFO_0000384")
# IBD
view("chr19:49206172_C_T","EFO_0003767")
# Systemic lupus erythematosus SLE
view("chr6:32424882_C_T","EFO_0002690")
# Rheumatoid arthritis
view("chr6:32424882_C_T","EFO_0000685")
# IgA nephropathy
view("chr6:32424882_C_T","EFO_0004194")
# Self-reported multiple sclerosis
view("chr6:32424882_C_T","EFO_0003885")
# Self-reported malabsorption or coeliac disease
view("chr6:32424882_C_T","EFO_0001060")
# Multiple sclerosis
view("chr6:32424882_C_T","EFO_0003885")
# Self-reported ankylosing spondylitis
view("chr6:32424882_C_T","EFO_0003898")
# Self-reported psoriasis
view("chr6:32424882_C_T","EFO_0000676")


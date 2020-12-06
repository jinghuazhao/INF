options(width=500)
library(dplyr)
library(pQTLtools)
catalogue <- "pQTL"; proxies <- "EUR"; p <- 5e-8; r2 <- 0.8; build <- 37
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE) %>% rename(INF1_rsid=rsid) %>% rename(Total=N)
INF1 <- within(left_join(metal,subset(inf1,select=-c(start,end))),{
                 hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                 HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)
               }) %>% rename(gene_gwas=gene) %>% rename(uniprot_gwas=uniprot)
r <- snpqueries(INF1[["INF1_rsid"]], catalogue="None", proxies="EUR", p=p, r2=r2, build=build)
INF1_aggr <- within(merge(INF1,with(r,snps)[c("snpid","hgnc")],by.x="MarkerName",by.y="snpid"),{
               gene_snpid <- paste0(hgnc,"-",MarkerName)
             }) %>% rename(gene=hgnc)
r <- snpqueries(INF1_aggr[["INF1_rsid"]], catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
ps <- subset(with(r,right_join(snps,results)),select=-c(ref_rsid, ref_hg19_coordinates, ref_hg38_coordinates,
                         ref_chr, ref_pos_hg19, ref_pos_hg38, ref_a1, ref_a2, ref_eur, ref_consequence,
                         ref_protein_position, ref_amino_acids, ref_ensembl, pos_hg19, pos_hg38, hg38_coordinates,
                         protein_position, amino_acids, ensembl, dprime, efo, n, n_studies, unit, direction))
f <- file.path(INF,"work","pQTL",paste0("INF1.merge.",catalogue))
save(INF1_aggr,r,ps,file=f)
ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
ips <- merge(within(INF1_aggr,{prot_snpid <- paste0(prot,"-",MarkerName)}),
             within(subset(ps,hgnc%in%INF1_aggr$gene),{gene_snpid <- paste0(hgnc,"-",snpid)}), by="gene_snpid",all.y=TRUE)
write.table(ips[c("uniprot_gwas","INF1_rsid","prot_snpid","cis.trans","proxy","r2","p","study","pmid","target.short","trait")],
      file="pQTL.log",quote=FALSE,row.names=FALSE,sep="\t")
simple <- ips%>%select(INF1_rsid,uniprot_gwas,prot_snpid,chr.x,chr.y,HLA,cis.trans,hgnc,proxy,r2,p,study,pmid,trait)
write.table(simple,file=file.path(INF,"work","pQTL","pQTL.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
repl <- with(SomaLogic160410,is.na(extGene))
SomaLogic160410[repl,"extGene"] <- SomaLogic160410[repl,"entGene"]
SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,extGene) %>% rename(hgnc=extGene)
pQTL <- dplyr::left_join(simple,SL)
INTERVAL <- subset(pQTL,pmid==29875488) %>%
            select(prot_snpid,chr.x,chr.y,chr,INF1_rsid,uniprot_gwas,HLA,cis.trans,
                   proxy,r2,p,study,pmid,chr,hgnc,trait,Target,TargetFullName)
write.table(INTERVAL,file=file.path(INF,"work","pQTL","pQTL-SomaLogic.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

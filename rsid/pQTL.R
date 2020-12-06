options(width=500)
library(dplyr)
library(pQTLtools)
catalogue <- "pQTL"; proxies <- "EUR"; p <- 5e-8; r2 <- 0.8; build <- 37
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
INF1 <- within(left_join(subset(metal,cis.trans=="trans"),subset(inf1,select=-c(start,end))),{
                 hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                 HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)
               }) %>% rename(INF1_rsid=rsid) %>% rename(Total=N) %>% rename(gene_gwas=gene) %>% rename(uniprot_gwas=uniprot)
r <- snpqueries(INF1[["INF1_rsid"]], catalogue="None", proxies="EUR", p=p, r2=r2, build=build)
snps <- with(r,snps)[c("snpid","hgnc")]
INF1_aggr <- within(merge(INF1,snps,by.x="MarkerName",by.y="snpid"), {gene_snpid <- paste0(hgnc,"-",MarkerName)}) %>% rename(gene=hgnc)
r <- snpqueries(INF1_aggr[["INF1_rsid"]], catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
ps <- with(r,right_join(snps,results))
f <- file.path(INF,"work",paste0("INF1.merge.",catalogue))
save(INF1_aggr,r,ps,file=f)
ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
ips <- subset(merge(INF1_aggr,within(subset(ps,hgnc%in%INF1_aggr$gene),{gene_snpid <- paste0(hgnc,"-",snpid)}),
                    by="gene_snpid",all.y=TRUE),select=-c(hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                        ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                        rsid, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                        dprime, efo, n, n_studies, unit, direction))
print(ips[c("prot","uniprot_gwas","INF1_rsid","gene_snpid","cis.trans","proxy","r2","study","pmid","target.short","trait")],row.names=FALSE,right=FALSE)
simple <- ips%>%select(INF1_rsid,prot,uniprot_gwas,target.short,gene_snpid,chr.x,chr.y,HLA,cis.trans,hgnc,proxy,r2,study,pmid,trait)
write.table(simple,file=file.path(INF,"work","pQTL.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
repl <- with(SomaLogic160410,is.na(extGene))
SomaLogic160410[repl,"extGene"] <- SomaLogic160410[repl,"entGene"]
SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,extGene) %>% rename(hgnc=extGene)
pQTL <- dplyr::left_join(simple,SL)
INTERVAL <- subset(pQTL,pmid==29875488) %>%
            select(gene_snpid,chr.x,chr.y,chr,INF1_rsid,prot,uniprot_gwas,HLA,cis.trans,
                   proxy,r2,study,pmid,chr,hgnc,trait,target.short,Target,TargetFullName)
write.table(INTERVAL,file=file.path(INF,"work","pQTL-SomaLogic.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

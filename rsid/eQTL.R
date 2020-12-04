options(width=500)

library(dplyr)
library(pQTLtools)

proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37
out <- "INF1.cis-GTEx.eQTL"

query <- function(rsid=INF1_aggr[["INF1_rsid"]])
for (catalogue in c("eQTL"))
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps,subset(results,study=="GTEx")))
  f <- file.path(INF,"work",out)
  save(INF1_aggr,r,ps,file=f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}
 
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
INF1 <- within(left_join(subset(metal,cis.trans=="cis"),inf1),{
                 hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                 HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)
               }) %>% rename(INF1_rsid=rsid) %>% rename(Total=N) %>% rename(gene_gwas=gene) %>% rename(uniprot_gwas=uniprot)
r <- snpqueries(INF1[["INF1_rsid"]], catalogue="None", proxies="EUR", p=p, r2=r2, build=build)
snps <- with(r,snps)[c("snpid","hgnc")] %>% rename(gene=hgnc)
# load(file.path(INF,"work","INF1.merge.trans.anno.rda"))
INF1_aggr <- within(merge(INF1,snps,by.x="MarkerName",by.y="snpid"), {gene_snpid <- paste0(gene,"-",MarkerName)})
# r <- snpqueries(snplist=with(trans,rsid),catalogue="None")
# m <- merge(trans,with(r,snps),by.x="MarkerName",by.y="snpid")
query()
f <- file.path(INF,"work",out)
load(f)
eQTL <- within(subset(ps,hgnc%in%INF1_aggr$gene), {gene_snpid <- paste0(hgnc,"-",snpid)}) %>%
              select(hgnc,ensembl,rsid,hg19_coordinates,a1,a2,eur,consequence,
              study,pmid,ancestry,year,tissue,exp_gene,exp_ensembl,beta,se,p,dataset,gene_snpid)
eQTL <- within(subset(eQTL,tissue!="Normal prepouch ileum"), {
  tissue <- gsub("Subcutaneous fat","Adipose subcutaneous",tissue)
  tissue <- gsub("Visceral abdominal fat","Adipose visceral omentum",tissue)
  tissue <- gsub("ba9","BA9",tissue)
  tissue <- gsub("ba24","BA24",tissue)
  tissue <- gsub("^Blood|Monocytes|Peripheral blood|Neutrophils|Peripheral blood monocytes|T cells|Whole Blood|Lymphoblastoid cell lines","Whole blood",tissue)
  tissue <- gsub("Breast tumors", "Breast mammary tissue", tissue)
  tissue <- gsub("Skin not sun exposed suprapubic|Skin sun exposed lower leg","Skin",tissue)
})
save(eQTL,file=file.path(INF,"work","cis-pQTL.eQTL.rda"))
eQTL_overlap <- subset(merge(INF1_aggr,eQTL,by="gene_snpid"),select=-c(Chromosome,Position))
write.table(eQTL_overlap,file=file.path(INF,"work","cis-pQTL_eQTL.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
tbl.cis <- with(within(subset(eQTL_overlap,cis.trans=="cis"),{rsidProts <- paste0(INF1_rsid," (",gene,")")}),table(rsidProts,tissue))
tbl.cis[tbl.cis>1] <- 1
sum(tbl.cis)
tbl.trans <- with(within(subset(eQTL_overlap,cis.trans=="trans"),{rsidProts <- paste0(INF1_rsid," (",gene,")")}),table(rsidProts,tissue))
tbl.trans[tbl.trans>1] <- 1
sum(tbl.trans)
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(INF1_rsid," (",gene,")")}),table(rsidProts,tissue))
tbl[tbl>1] <- 1
write.table(as.data.frame.matrix(tbl),file=file.path(INF,"work","cis-pQTL_eQTL_matrix.tsv"),quote=FALSE,sep="\t")
library(pheatmap)
pal <- colorRampPalette(c("white","red"))
col <- pal(3)
library(grid)
png(file.path(INF,"work","cis-pQTL_eQTL.png"),res=300,width=18,height=18,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(tbl, legend=FALSE, angle_col="45", color=col, width=8, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize=12)
setHook("grid.newpage", NULL, "replace")
grid.text("Tissue", y=-0.07, gp=gpar(fontsize=15))
grid.text("pQTL", x=-0.07, rot=90, gp=gpar(fontsize=15))
dev.off()

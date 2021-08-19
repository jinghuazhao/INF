query <- function(rsid=INF1_aggr[["INF1_rsid"]],catalogue="eQTL",keep=TRUE)
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps[c("snp","rsid","hg19_coordinates","a1","a2","consequence","hgnc","proxy","r2","snpid")],
                          subset(results,study=="GTEx",select=-c(ref_rsid, ref_hg19_coordinates, ref_hg38_coordinates, dprime,
                                                                 ref_a1, ref_a2, hg38_coordinates, efo, trait, probe, exp_ensembl,
                                                                 ancestry, year, direction,n,n_studies,unit,dataset))))
  f <- file.path(INF,"work",out)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  if (keep) save(INF1_aggr,r,ps,ips,file=f) else ps
}
 
options(width=500)
library(dplyr)
library(pQTLtools)
proxies <- "EUR"; p <- 5e-8; r2 <- 0.8; build <- 37; prefix <- "cis-pQTL-GTEx"; out <- paste0(prefix,".eQTL");
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
INF1 <- within(left_join(subset(metal,cis.trans=="cis"),subset(gap::inf1,select=-c(start,end))),{
                 hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                 HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)
               }) %>% rename(INF1_rsid=rsid, Total=N, gene_gwas=gene, uniprot_gwas=uniprot)
r <- snpqueries(INF1[["INF1_rsid"]], catalogue="None", proxies=proxies, p=p, r2=r2, build=build)
snps <- with(r,snps)[c("snpid","hgnc")] %>% rename(gene=hgnc)
INF1_aggr <- within(merge(INF1,snps,by.x="MarkerName",by.y="snpid"), {gene_snpid <- paste0(gene,"-",MarkerName)})
query()
f <- file.path(INF,"work",out)
load(f)
eQTL <- within(subset(ps,hgnc==exp_gene), {gene_snpid <- paste0(hgnc,"-",snpid)}) %>%
              select(hgnc,rsid,hg19_coordinates,a1,a2,proxy,r2,consequence,study,pmid,tissue,exp_gene,beta,se,p,gene_snpid)
eQTL <- within(subset(eQTL,tissue!="Normal prepouch ileum"), {
  tissue <- gsub("Subcutaneous fat","Adipose subcutaneous",tissue)
  tissue <- gsub("Visceral abdominal fat","Adipose visceral omentum",tissue)
  tissue <- gsub("ba9","BA9",tissue)
  tissue <- gsub("ba24","BA24",tissue)
  tissue <- gsub("^Blood|Monocytes|Peripheral blood|Neutrophils|Peripheral blood monocytes","Whole blood",tissue)
  tissue <- gsub("T cells|Whole Blood|Lymphoblastoid cell lines","Whole blood",tissue)
  tissue <- gsub("Breast tumors", "Breast mammary tissue", tissue)
  tissue <- gsub("Skin not sun exposed suprapubic|Skin sun exposed lower leg","Skin",tissue)
})
save(eQTL,file=file.path(INF,"work",paste0(prefix,".eQTL.rda")))
keep <- c("gene_snpid","INF1_rsid","prot", "gene", "uniprot_gwas", "gene_gwas", "cis.trans")
eQTL_overlap <- merge(INF1_aggr[keep],eQTL,by="gene_snpid")
write.table(eQTL_overlap,file=file.path(INF,"work",paste0(prefix,"_eQTL.tsv")),quote=FALSE,row.names=FALSE,sep="\t")
tbl.cis <- with(within(subset(eQTL_overlap,cis.trans=="cis"),{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),
                table(rsidProts,tissue))
tbl.cis[tbl.cis>1] <- 1
dim(tbl.cis)
sum(tbl.cis)
tbl.trans <- with(within(subset(eQTL_overlap,cis.trans=="trans"),{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),
                  table(rsidProts,tissue))
tbl.trans[tbl.trans>1] <- 1
dim(tbl.trans)
sum(tbl.trans)
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),table(rsidProts,tissue))
tbl[tbl>1] <- 1
write.table(as.data.frame.matrix(tbl),file=file.path(INF,"work",paste0(prefix,"_eQTL_matrix.tsv")),quote=FALSE,sep="\t")
library(pheatmap)
pal <- colorRampPalette(c("white","red"))
col <- pal(3)
library(grid)
png(file.path(INF,"work",paste0(prefix,"_eQTL.png")),res=300,width=18,height=18,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(tbl, legend=FALSE, angle_col="45", color=col, width=8, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize=22)
setHook("grid.newpage", NULL, "replace")
grid.text("Tissue", y=-0.07, gp=gpar(fontsize=15))
grid.text("pQTL", x=-0.07, rot=90, gp=gpar(fontsize=15))
dev.off()

chkList <- c("rs6827617", "rs4241577", "rs149278")
chkout <- query(chkList,keep=FALSE)
subset(chkout,hgnc==exp_gene,select=-c(hg19_coordinates,a1,a2,consequence,pmid))

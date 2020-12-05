query <- function(rsid=INF1_aggr[["INF1_rsid"]])
for (catalogue in c("eQTL","pQTL"))
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps,results))
  f <- paste0(file.path(INF,"work","INF1.merge."),catalogue)
  save(INF1_aggr,r,ps,file=f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}
 
options(width=500)
library(dplyr)
library(pQTLtools)
proxies <- "EUR"; p <- 5e-8; r2 <- 0.8; build <- 37
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
INF1 <- within(left_join(subset(metal,cis.trans=="trans"),inf1),{
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
f <- file.path(INF,"work","INF1.merge.eQTL")
load(f)
eQTL <- within(subset(ps,hgnc==exp_gene), {gene_snpid <- paste0(hgnc,"-",snpid)}) %>%
              select(hgnc,ensembl,rsid,hg19_coordinates,a1,a2,eur,consequence,
              study,pmid,ancestry,year,tissue,exp_gene,exp_ensembl,proxy,r2,beta,se,p,dataset,gene_snpid)
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
save(eQTL,file=file.path(INF,"work","pQTL.eQTL.rda"))
eQTL_overlap <- subset(merge(INF1_aggr,eQTL,by="gene_snpid"),select=-c(Chromosome,Position))
write.table(eQTL_overlap,file=file.path(INF,"work","pQTL_eQTL.tsv"),quote=FALSE,row.names=FALSE,sep="\t")
tbl.cis <- with(within(subset(eQTL_overlap,cis.trans=="cis"),{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),
                table(rsidProts,tissue))
tbl.cis[tbl.cis>1] <- 1
sum(tbl.cis)
tbl.trans <- with(within(subset(eQTL_overlap,cis.trans=="trans"),{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),
                  table(rsidProts,tissue))
tbl.trans[tbl.trans>1] <- 1
sum(tbl.trans)
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(INF1_rsid," (",gene_gwas,")")}),table(rsidProts,tissue))
tbl[tbl>1] <- 1
write.table(as.data.frame.matrix(tbl),file=file.path(INF,"work","pQTL_eQTL_matrix.tsv"),quote=FALSE,sep="\t")
library(pheatmap)
pal <- colorRampPalette(c("white","red"))
col <- pal(3)
library(grid)
png(file.path(INF,"work","pQTL_eQTL.png"),res=300,width=18,height=18,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(tbl, legend=FALSE, angle_col="45", color=col, width=8, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize=12)
setHook("grid.newpage", NULL, "replace")
grid.text("Tissue", y=-0.07, gp=gpar(fontsize=15))
grid.text("pQTL", x=-0.07, rot=90, gp=gpar(fontsize=15))
dev.off()
library(gap)
aux <- with(with(eQTL_overlap, cbind(inv_chr_pos_a1_a2(MarkerName)[c("chr","pos")],rsid,Allele1,Allele2,prot,HLA,cis.trans,tissue)), {
            flag <- (HLA==1)
            colId <- paste0(substr(chr,4,5),":",pos,"(",toupper(Allele1),"/",toupper(Allele2),")")
            colId[flag] <- paste0(colId[flag],"*")
            colLabel <- paste0(colId," (",prot,")")
            col <- rep("blue",nrow(eQTL_overlap))
            col[cis.trans=="cis"] <- "red"
            data.frame(colLabel,col,tissue)
          })
Col <- unique(aux[c("colLabel","col")])
rownames(Col) <- with(Col,colLabel)
library(gplots)
png(file.path(INF,"work","pQTL_eQTL_gplots.png"),height=35,width=40,units="cm",res=300)
heatmap.2(tbl, scale = "none", keysize=0.8, col = col, margin=c(20,20), trace = "none", key=FALSE,
          colCol=Col[colnames(tbl),"col"], dendrogram="none", density.info = "none", srtCol=45)
dev.off()

library(highcharter)
fntltp <- JS("function(){
 return this.series.xAxis.categories[this.point.x] + ' ' +
  this.series.yAxis.categories[this.point.y] + ':<br>' +
  Highcharts.numberFormat(this.point.value, 2);
}")
hc <- data.frame()
i <- 1
for(cn in colnames(tbl)) for(rn in rownames(tbl)) {
   hc[i,c("f1","f2","v")] <- c(cn,rn,tbl[rn,cn])
   i <- i + 1
}
n <- 4
stops <- data.frame(
  q = 0:n/n,
  c = c("#4287f5","grey","#ffffff","grey","#e32222"),
  stringsAsFactors = FALSE
  )
hc$f1 <- as.factor(hc$f1)
hc$f2 <- as.factor(hc$f2)
f1 <- levels(hc$f1)
highchart() %>%
  hc_title(text = "pQTLs and eQTL overlap",align="center")%>%
  hc_xAxis(categories = f1) %>%
  hc_yAxis(categories = hc$f2, reversed = TRUE)%>%
  hc_colorAxis(min = -1, max=1, stops=list_parse2(stops)) %>%
  hc_legend(align = "right",layout = "vertical",
            margin = 0,verticalAlign = "top",
            y = 30,symbolHeight = 200) %>%
  hc_tooltip(formatter = fntltp) %>%
  hc_add_series(data = hc, type = "heatmap",
                hcaes(x = f1,y = f2,value = v),
                dataLabels = list(enabled = FALSE))

rm(INF1_aggr,ps,r)
f <- file.path(INF,"work","INF1.merge.pQTL")
load(f)
ips <- subset(merge(INF1_aggr,within(subset(ps,hgnc%in%INF1_aggr$gene),{gene_snpid <- paste0(hgnc,"-",snpid)}),
                    by="gene_snpid",all.y=TRUE),select=-c(hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                        ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                        rsid, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                        dprime, efo, n, n_studies, unit, direction))
print(ips[c("prot","uniprot_gwas","INF1_rsid","gene_snpid","cis.trans","proxy","r2","study","pmid","target.short","trait")],row.names=FALSE,right=FALSE)
simple <- ips%>%select(INF1_rsid,prot,uniprot_gwas,target.short,gene_snpid,chr.x,chr.y,HLA,cis.trans,gene,proxy,r2,study,pmid,trait)
write.table(simple,file=file.path(INF,"work","pQTL.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# + INTERVAL SomaLogic data
repl <- with(SomaLogic160410,is.na(extGene))
SomaLogic160410[repl,"extGene"] <- SomaLogic160410[repl,"entGene"]
SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,extGene) %>% rename(gene=extGene)
pQTL <- dplyr::left_join(simple,SL)
INTERVAL <- subset(pQTL,pmid==29875488) %>%
            select(gene_snpid,chr.x,chr.y,chr,INF1_rsid,prot,uniprot_gwas,HLA,cis.trans,
                   proxy,r2,study,pmid,chr,gene,trait,target.short,Target,TargetFullName)
write.table(INTERVAL,file=file.path(INF,"work","pQTL-SomaLogic.tsv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# head -1 work/pQTL_eQTL_matrix.tsv | tr '\t' '\n' | grep -v signif
# ls *egene* | awk '{gsub(".v7.egenes.txt.gz","");print "#",NR" "$1}' | xsel -i

# 1 Adipose_Subcutaneous
# 2 Adipose_Visceral_Omentum
# 3 Adrenal_Gland
# 4 Artery_Aorta
# 5 Artery_Coronary
# 6 Artery_Tibial
# 7 Brain_Amygdala
# 8 Brain_Anterior_cingulate_cortex_BA24
# 9 Brain_Caudate_basal_ganglia
# 10 Brain_Cerebellar_Hemisphere
# 11 Brain_Cerebellum
# 12 Brain_Cortex
# 13 Brain_Frontal_Cortex_BA9
# 14 Brain_Hippocampus
# 15 Brain_Hypothalamus
# 16 Brain_Nucleus_accumbens_basal_ganglia
# 17 Brain_Putamen_basal_ganglia
# 18 Brain_Spinal_cord_cervical_c-1
# 19 Brain_Substantia_nigra
# 20 Breast_Mammary_Tissue
# 21 Cells_EBV-transformed_lymphocytes
# 22 Cells_Transformed_fibroblasts
# 23 Colon_Sigmoid
# 24 Colon_Transverse
# 25 Esophagus_Gastroesophageal_Junction
# 26 Esophagus_Mucosa
# 27 Esophagus_Muscularis
# 28 Heart_Atrial_Appendage
# 29 Heart_Left_Ventricle
# 30 Liver
# 31 Lung
# 32 Minor_Salivary_Gland
# 33 Muscle_Skeletal
# 34 Nerve_Tibial
# 35 Ovary
# 36 Pancreas
# 37 Pituitary
# 38 Prostate
# 39 Skin_Not_Sun_Exposed_Suprapubic
# 40 Skin_Sun_Exposed_Lower_leg
# 41 Small_Intestine_Terminal_Ileum
# 42 Spleen
# 43 Stomach
# 44 Testis
# 45 Thyroid
# 46 Uterus
# 47 Vagina
# 48 Whole_Blood

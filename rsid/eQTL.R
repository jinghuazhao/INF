query <- function(rsid=INF1_aggr[["INF1_rsid"]],catalogue="eQTL",keep=TRUE)
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps[c("snp","rsid","hg19_coordinates","a1","a2","consequence","hgnc","proxy","r2","snpid")],
                          subset(results,select=-c(ref_rsid, ref_hg19_coordinates, ref_hg38_coordinates, dprime,
                                                   ref_a1, ref_a2, hg38_coordinates, efo, trait, probe, exp_ensembl,
                                                   ancestry, year, direction,n,n_studies,unit,dataset))))
  f <- file.path(INF,"work",out)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  if (keep) save(INF1_aggr,r,ps,ips,file=f) else ps
}
 
options(width=500)
library(dplyr)
library(pQTLtools)
proxies <- "EUR"; p <- 5e-8; r2 <- 0.8; build <- 37; prefix <- "cis-pQTL"; out <- paste0(prefix,".eQTL");
INF <- Sys.getenv("INF")
metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE)
INF1 <- within(left_join(subset(metal,cis.trans=="cis"),subset(gap::inf1,select=-c(start,end))),{
                 hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                 HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)
               }) %>% rename(INF1_rsid=rsid) %>% rename(Total=N) %>% rename(gene_gwas=gene) %>% rename(uniprot_gwas=uniprot)
r <- snpqueries(INF1[["INF1_rsid"]], catalogue="None", proxies=proxies, p=p, r2=r2, build=build)
snps <- with(r,snps)[c("snpid","hgnc")] %>% rename(gene=hgnc)
# load(file.path(INF,"work","INF1.merge.trans.anno.rda"))
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
chkList <- c("rs6827617", "rs4241577", "rs149278")
chkout <- query(chkList,keep=FALSE)
subset(chkout,hgnc==exp_gene,select=-c(hg19_coordinates,a1,a2,consequence,pmid))
keep <- c("MarkerName","Allele1","Allele2","gene_snpid","INF1_rsid","prot", "gene", "uniprot_gwas", "gene_gwas", "cis.trans", "HLA")
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
cis_pQTL_eQTL <- eQTL_overlap %>% mutate(rsidProts=paste0(INF1_rsid," (",gene_gwas,")"),
                                          tissue_tags=paste0("P=",p,"[",tissue,",",study,",",pmid,"]")) %>%
                 group_by(rsidProts,tissue) %>%
                 summarize(rsidProts_tissue=paste(tissue_tags,collapse=";"))
cis_pQTL_eQTL_table <- with(cis_pQTL_eQTL,table(rsidProts,tissue))
for(row in rownames(cis_pQTL_eQTL_table)) for(col in colnames(cis_pQTL_eQTL_table))
  if (cis_pQTL_eQTL_table[row,col]==0) cis_pQTL_eQTL_table[row,col] <- "" else
  cis_pQTL_eQTL_table[row,col] <- subset(cis_pQTL_eQTL,rsidProts==row&tissue==col)[["rsidProts_tissue"]]
write.table(cis_pQTL_eQTL_table,file=file.path(INF,"work","cis_pQTL_eQTL_table.tsv"),sep="\t")

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
png(file.path(INF,"work",paste0(prefix,"_eQTL_gplots.png")),height=35,width=40,units="cm",res=300)
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

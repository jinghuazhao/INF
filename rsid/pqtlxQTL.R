options(width=200)

library(dplyr)
library(pQTLtools)

query <- function()
for (catalogue in c("eQTL","mQTL","pQTL"))
{
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  ps <- with(r,right_join(snps,results))
  f <- paste0(file.path(INF,"work","INF1.merge."),catalogue)
  save(INF1_aggr,r,ps,file=f)
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}
 
INF <- Sys.getenv("INF")
f <- paste0(file.path(INF,"work","INF1.merge."),"eQTL")
load(f)
INF1_aggr <- within(left_join(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),inf1),{
                      hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                      gene_snpid <- paste0(gene,"-",MarkerName)
                    }) %>% rename(INF1_rsid=rsid) %>% rename(Total=N)

rsid <- INF1_aggr[["INF1_rsid"]]
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37
# query()

eQTL <- within(ps,{gene_snpid <- paste0(hgnc,"-",snpid)}) %>% select(hgnc,ensembl,rsid,hg19_coordinates,a1,a2,eur,consequence,
                   study,pmid,ancestry,year,tissue,exp_gene,exp_ensembl,proxy,r2,beta,se,p,dataset,gene_snpid)
eQTL <- within(eQTL, {
  tissue <- gsub("ba9","BA9",tissue)
  tissue <- gsub("ba24","BA24",tissue)
  tissue <- gsub("^Blood|Monocytes|Peripheral blood|Neutrophils|Peripheral blood monocytes|T cells|Whole Blood","Whole blood",tissue)
})
save(eQTL,file="pQTL.eQTL.rda")
INF1_aggr <- within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)})
eQTL_overlap <- subset(merge(INF1_aggr,eQTL,by="gene_snpid"),select=-c(Chromosome,Position))
write.table(eQTL_overlap,file="pQTL_eQTL.tsv",quote=FALSE,row.names=FALSE,sep="\t")
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(INF1_rsid," (",prot,")")}),table(rsidProts,tissue))
tbl[tbl>1] <- 1
write.table(as.data.frame.matrix(tbl),file="pQTL_eQTL_matrix.tsv",quote=FALSE,sep="\t")
library(pheatmap)
pal <- colorRampPalette(c("white","red"))
col <- pal(3)
library(grid)
png("pQTL_eQTL.png",res=300,width=16,height=12,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(tbl, legend=FALSE, angle_col="45", color=col, width=8, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize=12)
setHook("grid.newpage", NULL, "replace")
grid.text("pQTL-protein", y=-0.07, gp=gpar(fontsize=15))
grid.text("Tissue", x=-0.07, rot=90, gp=gpar(fontsize=15))
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
png("pQTL_eQTL_gplots.png",height=35,width=40,units="cm",res=300)
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

ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates",all.y=TRUE),
              select=-c(snpid,
                        hg19_coordinates,hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                        ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                        rsid, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                        dprime, efo, n, n_studies, unit, direction))
SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,entGene) %>% rename(trait=TargetFullName)
pQTL <- dplyr::left_join(ips,SL)
write.table(pQTL,file=paste0(f,"-SomaLogic.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(subset(pQTL,UniProt%in%uniprot),file=paste0(f,"-ps.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

pQTL <- dplyr::left_join(subset(ips,pmid!=29875488),SL)
subset(pQTL[c("INF1_rsid","uniprot","UniProt","trait","ref_rsid","proxy","ref_a1","ref_a2","Allele1","Allele2","snp")],UniProt%in%uniprot)
# long and comprehensive list with less variables to fit screen
pQTL <- dplyr::left_join(subset(ips,pmid==29875488),SL)
subset(pQTL[c("INF1_rsid","uniprot","UniProt","trait","ref_rsid","proxy","ref_a1","ref_a2")],UniProt%in%uniprot)
gap::pvalue(-1.102/0.016)
gap::pvalue(-1.114/0.0157)
rs12075 <-c("P51671","P80162","P13500","P80075","P80098","Q99616")
subset(SomaLogic160410,UniProt%in%rs12075)
subset(SomaLogic160410,TargetFullName%in%c("Eotaxin","Corneodesmosin","C-C motif chemokine 14","C-C motif chemokine 26"))
subset(SomaLogic160410,TargetFullName=="SLAM family member 7")
subset(SomaLogic160410,UniProt=="P50591")
subset(SomaLogic160410,UniProt=="O14625")
subset(SomaLogic160410,UniProt=="P15692")

#genes <- scan("work/INF1.gene",what="")
#r <- genequeries(genes,catalogue="eQTL",build=37,p=5e-8,proxies="EUR",r2=0.8)

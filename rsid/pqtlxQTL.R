options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_aggr <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                      hg19_coordinates <- paste0("chr",Chromosome,":",Position)
                    }) %>% rename(INF1_rsid=rsid) %>% rename(Total=N) %>%
  select(Chromosome,Position,prot,uniprot,hg19_coordinates,MarkerName,INF1_rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
  group_by(Chromosome,Position,MarkerName,INF1_rsid,hg19_coordinates) %>% summarise(nprots <- n(),
           UniProts <- paste(uniprot,collapse=";"),
           prots <- paste(prot,collapse=";"),
           A1 <- paste(toupper(Allele1),collapse=";"),
           A2 <- paste(toupper(Allele2),collapse=";"),
           EAF <- paste(Freq1,collapse=";"),
           Effects <- paste(Effect,collapse=";"),
           SEs <- paste(StdErr,collapse=";"),
           log10P <- paste(log.P.,collapse=";"),
           cistrans <- paste(cis.trans,collapse=";")) %>% data.frame()

rsid <- INF1_aggr[["INF1_rsid"]]
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37

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
 
# pQTL <--> eQTL overlap
f <- paste0(file.path(INF,"work","INF1.merge."),"eQTL")
load(f)
eQTL <- ps %>% select(hgnc,ensembl,rsid,hg19_coordinates,a1,a2,eur,consequence,
                      study,pmid,ancestry,year,tissue,exp_gene,exp_ensembl,beta,se,p,dataset,snpid)
eQTL <- within(eQTL, {
  tissue <- gsub("ba9","BA9",tissue)
  tissue <- gsub("ba24","BA24",tissue)
  tissue <- gsub("^Blood|Monocytes|Peripheral blood|Neutrophils|Peripheral blood monocytes|T cells|Whole Blood","Whole blood",tissue)
})
save(eQTL,file="INF1.eQTL.rda")
INF1_aggr <- within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)})
eQTL_overlap <- subset(merge(INF1_aggr,eQTL,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
write.table(eQTL_overlap,file="INF1_eQTL.tsv",quote=FALSE,row.names=FALSE,sep="\t")
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(rsID," (",prots,")")}),table(rsidProts,tissue))
tbl[tbl>1] <- 1
library(pheatmap)
pal <- colorRampPalette(c("white","red"))
col <- pal(3)
## Create the heatmap:
library(grid)
png("INF1_eQTL.png",res=300,width=16,height=12,units="in")
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(tbl, legend=FALSE, angle_col="45", color=col, width=8, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize_row=6)
setHook("grid.newpage", NULL, "replace")
grid.text("pQTL-protein", y=-0.07, gp=gpar(fontsize=15))
grid.text("Tissue", x=-0.07, rot=90, gp=gpar(fontsize=15))
dev.off()
library(gap)
aux <- with(with(eQTL_overlap, cbind(inv_chr_pos_a1_a2(MarkerName)[c("chr","pos")],rsid,A1,A2,prots,HLA,cistrans,tissue)), {
            flag <- (HLA==1)
            colId <- paste0(substr(chr,4,5),":",pos,"(",A1,"/",A2,")")
            colId[flag] <- paste0(colId[flag],"*")
            colLabel <- paste0(colId," (",prots,")")
            col <- rep("blue",nrow(eQTL_overlap))
            col[cistrans=="cis"] <- "red"
            data.frame(colLabel,col,tissue)
          })
Col <- unique(aux[c("colLabel","col")])
rownames(Col) <- with(Col,colLabel)
library(gplots)
png("INF1_eQTL_gplots.png",height=35,width=40,units="cm",res=300)
heatmap.2(tbl, scale = "none", keysize=0.8, col = col, margin=c(20,20), trace = "none",
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

#genes <- scan("work/INF1.gene",what="")
#r <- genequeries(genes,catalogue="eQTL",build=37,p=5e-8,proxies="EUR",r2=0.8)
ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),
              select=-c(Chromosome, Position, EAF, Effects, SEs, nprots, snpid,
                        hg19_coordinates,hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                        ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                        rsid, chr, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                        dprime, efo, n, n_studies, unit, direction))
SL <- SomaLogic160410 %>% select(SOMAMER_ID,UniProt,Target,TargetFullName,chr,entGene) %>% rename(trait=TargetFullName)
pQTL <- dplyr::left_join(ips,SL)
write.table(pQTL,file=paste0(f,"-SomaLogic.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(subset(pQTL,UniProt%in%UniProts),file=paste0(f,"-ps.tsv"),row.names=FALSE,quote=FALSE,sep="\t")

pQTL <- dplyr::left_join(subset(ips,pmid!=29875488),SL)
subset(pQTL[c("rsID","UniProts","UniProt","trait","ref_rsid","proxy","ref_a1","ref_a2","A1","A2","snp")],UniProt%in%UniProts)
# long and comprehensive list with less variables to fit screen
pQTL <- dplyr::left_join(subset(ips,pmid==29875488),SL)
subset(pQTL[c("rsID","UniProts","UniProt","trait","ref_rsid","proxy","ref_a1","ref_a2")],UniProt%in%UniProts)
gap::pvalue(-1.102/0.016)
gap::pvalue(-1.114/0.0157)
rs12075 <-c("P51671","P80162","P13500","P80075","P80098","Q99616")
subset(SomaLogic160410,UniProt%in%rs12075)
subset(SomaLogic160410,TargetFullName%in%c("Eotaxin","Corneodesmosin","C-C motif chemokine 14","C-C motif chemokine 26"))
subset(SomaLogic160410,TargetFullName=="SLAM family member 7")
subset(SomaLogic160410,UniProt=="P50591")
subset(SomaLogic160410,UniProt=="O14625")
subset(SomaLogic160410,UniProt=="P15692")

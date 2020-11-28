options(width=200)

library(dplyr)
library(pQTLtools)

INF <- Sys.getenv("INF")
INF1_aggr <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                     hg19_coordinates=paste0("chr",Chromosome,":",Position)}) %>% rename(INF1_rsid=rsid) %>% rename(Total=N) %>%
  select(Chromosome,Position,prot,uniprot,hg19_coordinates,MarkerName,INF1_rsid,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans) %>%
  group_by(Chromosome,Position,MarkerName,INF1_rsid,hg19_coordinates) %>% summarise(nprots=n(),
           UniProts=paste(uniprot,collapse=";"),
           prots=paste(prot,collapse=";"),
           A1=paste(toupper(Allele1),collapse=";"),
           A2=paste(toupper(Allele2),collapse=";"),
           EAF=paste(Freq1,collapse=";"),
           Effects=paste(Effect,collapse=";"),
           SEs=paste(StdErr,collapse=";"),
           log10P=paste(log.P.,collapse=";"),
           cistrans=paste(cis.trans,collapse=";")) %>% data.frame()

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
  ips <- subset(merge(INF1_aggr,ps,by="hg19_coordinates"),
                select=-c(Chromosome, Position, EAF, Effects, SEs, nprots, snpid,
                          hg19_coordinates,hg38_coordinates,ref_hg19_coordinates,ref_hg38_coordinates,
                          ref_pos_hg19, ref_pos_hg38, ref_protein_position, ref_amino_acids, ref_ensembl,
                          rsid, chr, pos_hg19, pos_hg38, protein_position, amino_acids, ensembl,
                          dprime, efo, n, n_studies, unit, direction))
  write.table(ips,file=paste0(f,".tsv"),row.names=FALSE,quote=FALSE,sep="\t")
}

# pQTL <--> eQTL overlap
f <- paste0(file.path(INF,"work","INF1.merge."),"eQTL")
eQTL <- ps %>%
        select(hgnc,ensembl,rsid,hg19_coordinates,a1,a2,eur,consequence,
               ensembl,hgnc,study,pmid,ancestry,year,tissue,exp_gene,exp_ensembl,beta,se,p,dataset,snpid)
eQTL <- within(eQTL, {
  tissue <- gsub("ba9","BA9",tissue)
  tissue <- gsub("ba24","BA24",tissue)
  tissue <- gsub("^Blood|Monocytes|Peripheral blood|Neutrophils|Peripheral blood monocytes|T cells|Whole Blood","Whole blood",tissue)
})
save(eQTL,file="INF1.eQTL.rda")
INF1_aggr <- within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)})
eQTL_overlap <- subset(merge(INF1_aggr,eQTL,by="hg19_coordinates"),select=-c(hg19_coordinates,Chromosome,Position))
tbl <- with(within(eQTL_overlap,{rsidProts <- paste0(rsID," (",prots,")")}),table(tissue,rsidProts))
write.table(eQTL_overlap,file="INF1_eQTL.tsv",quote=FALSE,row.names=FALSE,sep="\t")

library(pheatmap)
pal <- colorRampPalette(c("white","blue", "red"))
col <- pal(20)
pheatmap(tbl,
         color = col,
         legend = TRUE,
         main = "Olink pQTLs overlapping with eQTLs across tissues",
         angle_col = "45",
         filename = "INF1_eQTL.png",
         width = 16,
         height = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellheight = 20,
         cellwidth = 20,
         fontsize = 11)
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

#genes <- scan("work/INF1.gene",what="")
#r <- genequeries(genes,catalogue="eQTL",build=37,p=5e-8,proxies="EUR",r2=0.8)

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

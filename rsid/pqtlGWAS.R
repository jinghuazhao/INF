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
p <- 5e-8
r2 <- 0.8
build <- 37

r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
lapply(r,dim)
snps_results <- with(r,right_join(snps,results))
ps <- subset(snps_results,select=-c(hg38_coordinates,ref_hg38_coordinates,pos_hg38,ref_pos_hg38,dprime))
sarcoidosis <- with(ps,grep("sarcoidosis",trait))
ps[sarcoidosis,"efo"] <- "Orphanet_797"
save(INF1_aggr,r,ps,file=file.path(INF,"work","INF1.merge.GWAS"))

metal <- subset(within(INF1_metal,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
                select=-c(Chromosome,Position,INF1_rsid,Direction))
aggr <- subset(within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
               select=-c(Chromosome,Position,INF1_rsid))
immune_infection <- read.delim("doc/immune.efos.txt",as.is=TRUE)
immune_infection_efo <- with(immune_infection,gsub(":","_",id))
long <- cbind(merge(metal,subset(ps,efo%in%immune_infection_efo),by="hg19_coordinates"),infection=0)
short <- cbind(merge(aggr,subset(ps,efo%in%immune_infection_efo),by="hg19_coordinates"),infection=0)
infection_efo <- with(subset(immune_infection,infection==1),gsub(":","_",id))
long[with(long,efo)%in%infection_efo,"infection"] <- 1
short[with(short,efo)%in%infection_efo,"infection"] <- 1

require(openxlsx)
xlsx <- "work/pqtl-immune_infection.xlsx"
wb <- createWorkbook(xlsx)
addWorksheet(wb, "METAL")
writeDataTable(wb, "METAL", subset(INF1_metal,select=-c(INF1_rsid,hg19_coordinates)))
addWorksheet(wb, "ps")
writeDataTable(wb, "ps", ps)
addWorksheet(wb, "EFO")
writeDataTable(wb, "EFO", immune_infection)
addWorksheet(wb, "long")
writeDataTable(wb, "long", 
               subset(long,select=-c(hg19_coordinates,ref_hg19_coordinates,ref_protein_position,ref_amino_acids,
                      snp,snpid,protein_position,amino_acids)))
addWorksheet(wb, "short")
writeDataTable(wb, "short", 
               subset(short,select=-c(hg19_coordinates,ref_hg19_coordinates,ref_protein_position,ref_amino_acids,
                      snp,snpid,protein_position,amino_acids)))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

view <- function(id,efoid,
                 v=c("MarkerName","Allele1","Allele2","a1","a2","efo","ref_a1","ref_a2","proxy","r2",
                     "beta","se","p","trait","ancestry","pmid","study"))
{
  options(width=200)
  cat(id,efoid,"\n")
  d <- subset(short[v],MarkerName==id & efo==efoid)
  subset(d,select=-c(MarkerName,efo))
}
# Rheumatoid arthritis
view("chr1:154426970_A_C","EFO_0000685")
# T1D
view("chr12:111884608_C_T", "EFO_0001359")
# Primary sclerosing cholangitis
view("chr12:111884608_C_T", "EFO_0004268")
# celiac
view("chr12:111884608_C_T", "EFO_0001060")
# Hypothyroidism
view("chr12:111884608_C_T", "EFO_0004705")
# Inflammatory bowel disease
view("chr12:111884608_C_T", "EFO_0003767")
# Primary biliary cirrhosis
view("chr12:111884608_C_T", "EFO_1001486")
# Allergic disease asthma hay fever or eczema
view("chr12:111932800_C_T", "EFO_0003785")
# Celiac
view("chr12:112007756_C_T","EFO_0001060")
# Crohn's
view("chr19:49206172_C_T","EFO_0000384")
# IBD
view("chr19:49206172_C_T","EFO_0003767")
# Rheumatoid arthritis
view("chr6:32424882_C_T","EFO_0000685")
# Self-reported ankylosing spondylitis
view("chr6:32424882_C_T","EFO_0003898")
# Multiple sclerosis
view("chr6:32424882_C_T","EFO_0003885")
# Self-reported psoriasis
view("chr6:32424882_C_T","EFO_0000676")
# Systemic lupus erythematosus SLE
view("chr6:32424882_C_T","EFO_0002690")
# Self-reported malabsorption or coeliac disease
view("chr6:32424882_C_T","EFO_0001060")
# IgA nephropathy
view("chr6:32424882_C_T","EFO_0004194")
# Self-reported sarcoidosis
view("chr6:32424882_C_T","Orphanet_797")
# Primary sclerosing cholangitis
view("chr6:32424882_C_T","EFO_0004268")

xlsx <- "work/pqtl-immune_infection_edited.xlsx"
pqtl_immune_infection <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:51), rows=c(1:220))
v=c("prots","MarkerName","Effects","Allele1","Allele2","rsid","a1","a2","efo","ref_rsid","ref_a1","ref_a2","proxy","r2",
    "HLA","infection","beta","se","p","trait","n_cases","n_controls","unit","ancestry","pmid","study","Keep","Switch")
mat <- within(subset(pqtl_immune_infection,infection==0 & Keep==1)[v],
{
  flag <- (HLA==1)
  prefix <- rsid
  prefix[flag] <- paste0(rsid[flag],"*")
  rsidProts <- paste0(prefix," (",prots,")")
  trait_shown <- gsub("Self-reported |Other |Doctor diagnosed ","",trait)
  trait_shown <- gsub("asthma |Allergic disease asthma hay fever or eczema","Allergic disease",trait_shown)
  trait_shown <- gsub("celiac disease|Celiac disease","malasorption or celiac disease",trait_shown)
  trait_shown <- gsub("malabsorption or coeliac disease","malasorption or celiac disease",trait_shown)
  trait_shown <- gsub("systemic lupus erythematosis or sle|Systemic lupus erythematosus SLE","systemic lupus erythematosus",trait_shown)
  trait_shown <- gsub("\\b(^[a-z])","\\U\\1",trait_shown, perl=TRUE)
  positiveEffects <- sign(as.numeric(Effects))
  positiveEffects[is.na(as.numeric(Effects))] <- 1
  qtl_direction <- sign(as.numeric(beta))
  qtl_direction[unit=="-"] <- 0.5*(runif(1)>0)
  qtl_direction[positiveEffects==-1] <- -qtl_direction[positiveEffects==-1]
  qtl_direction[!is.na(Switch)] <- -qtl_direction[!is.na(Switch)]
  efoTraits <- paste0(gsub("_",":",efo)," (",trait_shown,")")
})
# all beta's are NAs when unit=="-"
subset(mat[c("study","pmid","unit","beta","Keep")],unit=="-")
# all studies with risk difference were UKBB
subset(mat[c("study","pmid","unit","beta","n_cases","n_controls","Keep")],unit=="risk diff")

rxc <- with(mat,table(efoTraits,rsidProts))
indices <- mat[c("efoTraits","rsidProts","qtl_direction")]
# rxc_indices <- rxc[indices[,1],indices[,2]]
# ord <- order(indices[,1],indices[,2])
# rxc_indices[rxc_indices==1] <- indices[ord,3]
for(cn in colnames(rxc)) for(rn in rownames(rxc)) {
   s <- subset(indices,efoTraits==rn & rsidProts==cn)
   qd <- s[["qtl_direction"]]
   if(length(qd)>1) stop("duplicates")
   class(qd) <- "numeric"
   if(nrow(s)>0 & !is.na(qd[1])) rxc[rn,cn] <- qd[1]
}

library(pheatmap)
col <- colorRampPalette(c("#4287f5","grey","#ffffff","grey","#e32222"))(5)
pheatmap(rxc,
         color = col,
         legend = TRUE,
         main = "Olink pQTLs overlapping with QTLs for immune outcomes",
         angle_col = "45",
         filename = "INF1_pQTL_immune_qtl_unclustered.png",
         width = 16,
         height = 10,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellheight = 20,
         cellwidth = 20)

pheatmap(rxc,
         color = col,
         legend = TRUE,
         main = "Olink pQTLs overlapping with QTLs for immune outcomes",
         angle_col = "45",
         filename = "INF1_pQTL_immune_qtl.png",
         width = 16,
         height = 10,
         treeheight_row = 100,
         treeheigh_col = 100,
         cellheight = 20,
         cellwidth = 20)

obsolete <- function()
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

# A test of colorRampPalette
  YlOrBr <- c("#4287f5","grey","#ffffff","grey","#e32222")
  filled.contour(volcano,color.palette = colorRampPalette(YlOrBr, space = "Lab"), asp = 1)
# Colouring for the dendrogram
  library(dendextend)
  Rowv <- rxc %>% scale %>% dist %>% hclust %>% as.dendrogram %>%
     set("branches_k_color", k = 3) %>% set("branches_lwd", 1.2) %>%
     ladderize
  Colv <- rxc %>% scale %>% t %>% dist %>% hclust %>% as.dendrogram %>%
     set("branches_k_color", k = 2, value = c("orange", "blue")) %>%
     set("branches_lwd", 1.2) %>%
     ladderize
# stats
  heatmap(scale(rxc), scale = "none")
  heatmap(scale(rxc), Rowv = Rowv, Colv = Colv, scale = "none")
  library("gplots")
  heatmap.2(rxc, scale = "none", col = colorpanel(3, "blue", "white", "red"), trace = "none", density.info = "none")
  heatmap.2(scale(rxc), scale = "none", col = bluered(100), Rowv = Rowv, Colv = Colv, trace = "none", density.info = "none")
}

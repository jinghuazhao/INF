options(width=200)

suppressMessages(library(dplyr))
suppressMessages(library(gap))
suppressMessages(library(pQTLtools))

inf1_prot <- vector()
for(i in 1:92) inf1_prot[inf1[i,"prot"]] <- mutate(inf1[i,],target.short=if_else(!is.na(alt_name),alt_name,target.short))[["target.short"]]
INF <- Sys.getenv("INF")
INF1_metal <- within(read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE),{
                    hg19_coordinates=paste0("chr",Chromosome,":",Position)}) %>%
                    rename(INF1_rsid=rsid, Total=N) %>%
                    left_join(pQTLdata::inf1[c("prot","gene","target.short","alt_name")]) %>%
                    mutate(target.short=if_else(!is.na(alt_name),alt_name,target.short)) %>%
                    select(-alt_name)
INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,target.short,gene,hg19_coordinates,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans,INF1_rsid) %>%
  group_by(Chromosome,Position,MarkerName,INF1_rsid,hg19_coordinates) %>%
  summarise(nprots=n(),
            prots=paste(target.short,collapse=";"),
            Allele1=paste(toupper(Allele1),collapse=";"),
            Allele2=paste(toupper(Allele2),collapse=";"),
            EAF=paste(Freq1,collapse=";"),
            Effects=paste(Effect,collapse=";"),
            SEs=paste(StdErr,collapse=";"),
            log10P=paste(log.P.,collapse=";"),
            cistrans=paste(cis.trans,collapse=";")) %>% data.frame()
INF1_aggr_save <- INF1_aggr
setwd(INF)
rsid <- INF1_aggr[["INF1_rsid"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37
efo_diseases <- read.table(file.path(INF,"ebi","efo-3.26.0","efo_diseases.csv"),col.names=c("efo","disease"),as.is=TRUE,sep=",") %>%
                mutate(efo=gsub(":", "_", efo))

if (FALSE) {
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  snps_results <- with(r,right_join(snps,results))
  ps <- subset(snps_results,select=-c(hg38_coordinates,ref_hg38_coordinates,pos_hg38,ref_pos_hg38,dprime))
  sarcoidosis <- with(ps,grep("sarcoidosis",trait))
  ps[sarcoidosis,"efo"] <- "Orphanet_797"
  INF1_aggr <- INF1_aggr_save
  save(INF1_aggr,r,ps,file=file.path(INF,"work","INF1.merge.GWAS"))
} else {
  load(file.path(INF,"work","INF1.merge.GWAS"))
  INF1_aggr <- INF1_aggr_save
}
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
xlsx <- file.path(INF,"work","pqtl-immune_infection.xlsx")
if (FALSE)
{
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
}

view <- function(id,efoid,
                 v=c("MarkerName","Allele1","Allele2","a1","a2","efo","ref_a1","ref_a2","proxy","r2",
                     "beta","se","p","trait","ancestry","pmid","study"))
{
  options(width=200)
  cat(id,efoid,"\n")
  d <- subset(short[v],MarkerName==id & efo==efoid)
  subset(d,select=-c(MarkerName,efo))
}

imd <- function()
{
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
  short <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:51), rows=c(1:220))
  for(i in 1:nrow(short))
  {
    nprots <- short[i,"nprots"]
    ij <- unlist(strsplit(short[i,"prots"],";"))
    for(j in 1:nprots) short[i, "prots"] <- gsub(ij[j],inf1_prot[ij[j]],short[i,"prots"])
  }
  v <- c("prots","hgnc","MarkerName","cistrans","Effects","Allele1","Allele2","rsid","a1","a2","efo",
         "ref_rsid","ref_a1","ref_a2","proxy","r2",
         "HLA","infection","beta","se","p","trait","n_cases","n_controls","unit","ancestry","pmid","study","Keep","Switch")
  mat <- within(subset(short,infection==0 & Keep==1)[v],
  {
    flag <- (HLA==1)
    prefix <- paste0(prots,"-",rsid)
    prefix[flag] <- paste0(prefix[flag],"*")
    rsidProts <- paste0(stringr::str_pad(gsub("chr|:[0-9]*|_[A-Z]*","",MarkerName), width=2, side="left", pad="0"),"-", prefix," (",hgnc,")")
    trait_shown <- gsub("Self-reported |Other |Doctor diagnosed ","",trait)
    trait_shown <- gsub("asthma |Allergic disease asthma hay fever or eczema","Allergic disease",trait_shown)
    trait_shown <- gsub("celiac disease|Celiac disease","malasorption or celiac disease",trait_shown)
    trait_shown <- gsub("malabsorption or coeliac disease","malasorption or celiac disease",trait_shown)
    trait_shown <- gsub("systemic lupus erythematosis or sle|Systemic lupus erythematosus SLE","systemic lupus erythematosus",trait_shown)
    trait_shown <- gsub("\\b(^[a-z])","\\U\\1",trait_shown, perl=TRUE)
  # positiveEffects <- sign(as.numeric(Effects))
  # it happened that all NA's have beta>0 from multiple proteins
  # positiveEffects[is.na(as.numeric(Effects))] <- 1
    qtl_direction <- sign(as.numeric(beta))
  # we have >1 with beta=NA, so settle on either side of zero
  # qtl_direction[unit=="-"] <- 0.5*(runif(1)>0)
  # qtl_direction[positiveEffects==-1] <- -qtl_direction[positiveEffects==-1]
    qtl_direction[!is.na(Switch)] <- -qtl_direction[!is.na(Switch)]
  # efoTraits <- paste0(gsub("_",":",efo)," (",trait_shown,")")
    efoTraits <- paste0(trait_shown)
  })
  rxc <- with(mat,table(efoTraits,rsidProts))
  indices <- mat[c("efoTraits","rsidProts","qtl_direction")]
  if (FALSE) {
    add_entry <- data.frame(efoTraits="Multiple sclerosis",rsidProts="12-rs2364485 [TNFB](LTA)",qtl_direction=1)
    indices_new <- rbind(indices,add_entry)
    rxc <- with(indices_new,table(efoTraits,rsidProts))
  }
  for(cn in colnames(rxc)) for(rn in rownames(rxc)) {
     cnrn <- subset(indices,efoTraits==rn & rsidProts==cn)
     if(nrow(cnrn)==0) next
     qd <- as.numeric(cnrn[["qtl_direction"]])
     if(length(qd)>1) stop("duplicates")
     rxc[rn,cn] <- qd[1]
  }
  rxc[is.na(rxc)] <- 0
  rxc
}

rxc <- imd()

SF <- function(rxc, f="SF-pQTL-IMD-GWAS.png", ch=21, cw=21, h=13, w=18, ylab="Immune-mediated outcomes")
{
  library(pheatmap)
  col <- colorRampPalette(c("#4287f5","#ffffff","#e32222"))(3)
  library(grid)
  png(file.path(INF,f),res=300,width=w,height=h,units="in")
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  colnames(rxc) <- gsub("^[0-9]*-","",colnames(rxc))
  pheatmap(rxc, legend=FALSE, angle_col="315", border_color="black", color=col, cellheight=ch, cellwidth=cw, cluster_rows=TRUE, cluster_cols=FALSE, fontsize=16)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Protein(s)-pQTL (gene)", y=-0.07, gp=gpar(fontsize=15))
  grid.text(ylab, x=-0.07, rot=90, gp=gpar(fontsize=15))
  dev.off()
}

# INF1_pQTL_immune_qtl_unclustered.png
SF(rxc)

gwas <- function()
{
  short <- merge(aggr,ps,by="hg19_coordinates") %>%
           filter(efo %in% pull(efo_diseases,efo)) %>%
           left_join(efo_diseases)
  v <- c("prots","hgnc","MarkerName","cistrans","Effects","Allele1","Allele2","rsid","a1","a2","efo",
         "ref_rsid","ref_a1","ref_a2","proxy","r2",
         "HLA","beta","se","p","disease","n_cases","n_controls","unit","ancestry","pmid","study")
  mat <- within(short[v],
  {
    flag <- (HLA==1)
    prefix <- paste0(prots,"-",rsid)
    prefix[flag] <- paste0(prefix[flag],"*")
    rsidProts <- paste0(prefix," (",hgnc,")")
    efoTraits <- gsub("\\b(^[a-z])","\\U\\1",disease,perl=TRUE)
    qtl_direction <- sign(as.numeric(beta))
  })
  combined <- group_by(mat,efoTraits,rsidProts,desc(n_cases)) %>%
              summarize(direction=paste(qtl_direction,collapse=";"),
                        betas=paste(beta,collapse=";"),
                        units=paste(unit,collapse=";"),
                        studies=paste(study,collapse=";"),
                        diseases=paste(disease,collapse=";"),
                        cases=paste(n_cases,collapse=";")
                       ) %>%
              data.frame()
  rxc <- with(combined,table(efoTraits,rsidProts))
  for(cn in colnames(rxc)) for(rn in rownames(rxc)) {
     cnrn <- subset(combined,efoTraits==rn & rsidProts==cn)
     if(nrow(cnrn)==0) next
     rxc[rn,cn] <- as.numeric(unlist(strsplit(cnrn[["direction"]],";"))[1])
  }
  rxc[is.na(rxc)] <- 0
  # all beta's are NAs when unit=="-"
  subset(mat[c("study","pmid","unit","beta","qtl_direction")],unit=="-")
  # all studies with risk difference were UKBB
  subset(mat[c("study","pmid","unit","beta","n_cases","n_controls","qtl_direction")],unit=="risk diff")
  write.table(mat,file=file.path(INF,"work","pQTL-disease-GWAS.csv"),row.names=FALSE,quote=FALSE,sep=",")
  write.table(combined,file=file.path(INF,"work","pQTL-disease-GWAS-combined.csv"),row.names=FALSE,quote=FALSE,sep=",")
  rxc
}

rxc <- gwas()

# GWAS traits and diseases
SF(rxc,f="SF-pQTL-GWAS.png",ch=21,cw=21,h=30,w=20,ylab="GWAS diseases")

obsolete <- function()
{
  pheatmap(rxc,
           color = col,
           legend = TRUE,
           main = "",
           angle_col = "315",
           filename = "INF1_pQTL_immune_qtl.png",
           width = 17,
           height = 11,
           treeheight_row = 100,
           treeheigh_col = 100,
           cellheight = 20,
           cellwidth = 20,
           fontsize_row = 14,
           fontsize = 13)

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
    hc_title(text = "pQTLs and Immune-related Diseases",align="center")%>%
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

  library(gap)
  aux <- with(with(mat, cbind(inv_chr_pos_a1_a2(MarkerName)[c("chr","pos")],rsid,Allele1,Allele2,prots,HLA,cistrans,efoTraits,qtl_direction)), {
              flag <- (HLA==1)
  # a bit too specific here as they involve too many proteins nevertheless each combination with the same effect allele
              Allele1[8:13] <- "T"
              Allele2[8:13] <- "C"
              Allele1[25] <- "C"
              Allele2[25] <- "G"
              colId <- paste0(substr(chr,4,5),":",pos,"(",Allele1,"/",Allele2,")")
              colId[flag] <- paste0(colId[flag],"*")
              colLabel <- paste0(colId," (",prots,")")
              col <- rep("blue",nrow(mat))
              col[cistrans=="cis"] <- "red"
              data.frame(colLabel,col,efoTraits,qtl_direction)
            })
  Col <- unique(aux[c("colLabel","col")])
  rownames(Col) <- with(Col,colLabel)

  RXC <- with(aux,table(efoTraits,colLabel))
  indices <- aux[c("efoTraits","colLabel","qtl_direction")]
  for(cn in colnames(RXC)) for(rn in rownames(RXC)) {
     s <- subset(indices,efoTraits==rn & colLabel==cn)
     qd <- s[["qtl_direction"]]
     if(length(qd)>1) stop("duplicates")
     class(qd) <- "numeric"
     if(nrow(s)>0 & !is.na(qd[1])) RXC[rn,cn] <- qd[1]
  }

  library(gplots)
  png("INF1_pQTL_immune_gplots.png",height=35,width=40,units="cm",res=300)
  heatmap.2(RXC, scale = "none", keysize=0.8, col = colorpanel(5, "blue", "white", "red"), margin=c(20,20), trace = "none",
            colCol=Col[colnames(RXC),"col"], dendrogram="none", density.info = "none", srtCol=45)
  dev.off()

  efo_list_immune <- subset(read.csv("work/efo_list_annotated.csv",as.is=TRUE),immune_mediated==1)
  isd1 <- merge(aggr,subset(ps,efo%in%with(efo_list_immune,EFO)),by="hg19_coordinates")
  write.table(isd1,file="isd1.tsv",row.names=FALSE,quote=FALSE,sep="\t")

  load(file.path(INF,"files","efo.rda"))
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
# gplots
  heatmap.2(scale(rxc), scale = "none", col = bluered(100), Rowv = Rowv, Colv = Colv, trace = "none", density.info = "none")
  options(width=200)
  IL_12Bcis <- ieugwasr::phewas("rs10076557")
  data.frame(IL_12Bcis)
}

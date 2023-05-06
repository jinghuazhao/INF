options(width=2000)
require(openxlsx)
url <- "https://jhz22.user.srcf.net/INF1.latest.xlsx"
INF <- Sys.getenv("INF")
url <- file.path(INF,"work","INF1.latest.xlsx")
read.sheet <- function(sheet,cols,rows) read.xlsx(url,sheet=sheet,colNames=TRUE,cols=cols,rows=rows,skipEmptyRows=TRUE)
suppressMessages(library(dplyr))
suppressMessages(library(gap))
require(stringr)
gap_inf1 <- gap.datasets::inf1[c("uniprot", "prot", "target.short")] %>%
            mutate(target.short=gsub("MCP-1","CCL2",target.short)) %>%
            mutate(target.short=gsub("MCP-2","CCL8",target.short)) %>%
            mutate(target.short=gsub("MCP-3","CCL7",target.short)) %>%
            mutate(target.short=gsub("MCP-4","CCL13",target.short))
 summary <- read.sheet("Summary", 1:2, 2:36)
    inf1 <- subset(read.sheet("INF1", 1:12, 2:94),uniprot!="P23560") %>%
            select(-panel) %>%
            setNames(c("Target","Protein","UniProt","onMultiplePanel","onPanels","hgncSymbol","Chromosome","Start","End","olinkId","altUniProt")) %>%
            mutate(Protein=gsub("MCP-1","CCL2",Protein),
                   Protein=gsub("MCP-2","CCL8",Protein),
                   Protein=gsub("MCP-3","CCL7",Protein),
                   Protein=gsub("MCP-4","CCL13",Protein)) %>%
            select(UniProt,Target,Protein,onMultiplePanel,onPanels,hgncSymbol,Chromosome,Start,End,olinkId,altUniProt)
 studies <- read.sheet("Studies", 1:3, 2:14) %>%
            filter(Name!="MadCam")
 studies <- bind_rows(studies,data.frame(Name="Total",Design="",Size=sum(studies$Size))) %>%
            mutate(Size=formatC(Size,format="f",big.mark=",",digits=0,width=11))
   pqtls <- merge(read.sheet("pQTLs", 1:21, 2:182),gap_inf1[c("prot","target.short")],by="prot")
interval <- merge(within(read.sheet("INTERVAL", 1:12, 2:29),{Protein <- gsub(" ", "", Protein)}),
                  gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short,r2=as.character(r2),p=as.character(p),PMID=as.character(PMID)) %>% select(-target.short)
      os <- merge(read.sheet("OtherStudies", 1:12, 2:102),gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short,r2=as.character(r2),p=as.character(p),PMID=as.character(PMID)) %>% select(-target.short)
    cvd1 <- merge(read.sheet("CVD1", 1:12, 2:61),gap_inf1[c("prot","target.short")],by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short,r2=as.character(r2),p=as.character(p),PMID=as.character(PMID)) %>% select(-target.short)
 fenland <- read.delim(file.path(INF,"Fenland","Fenland.tsv")) %>%
            mutate(r2=as.character(r2),p=as.character(p),PMID=as.character(PMID),Comment=as.character(Comment))
  decode <- read.delim(file.path(INF,"deCODE","deCODE.tsv")) %>%
            mutate(r2=as.character(r2),p=as.character(p),PMID=as.character(PMID),Comment=as.character(Comment))
    aric <- read.delim(file.path(INF,"ARIC","ARIC.tsv")) %>%
            mutate(r2=as.character(r2),p=as.character(p),PMID=as.character(PMID),Comment=as.character(Comment))
    ages <- read.delim(file.path(INF,"AGES","AGES.tsv")) %>%
            mutate(r2=as.character(r2),p=as.character(p),PMID=as.character(PMID),Comment=as.character(Comment))
     ngs <- read.delim(file.path(INF,"ukb","NGS.tsv")) %>%
            mutate(r2=as.character(r2),p=as.character(p),PMID=as.character(PMID),Comment=as.character(Comment))
aristotl <- merge(read.sheet("ARISTOTLE", 1:14, 2:182), gap_inf1[c("prot","target.short")], by.x="Protein", by.y="prot") %>%
            filter(!is.na(POS)) %>%
            rename(rsid=SNPID,Chromosome=CHR,Position=POS,Allele1=EFFECT_ALLELE,Allele2=REFERENCE_ALLELE,
                   EAF=CODE_ALL_FQ,b=BETA,P=PVAL,Info=RSQ_IMP,Imputed=IMP) %>%
            mutate(Protein=target.short,Position=formatC(Position,format="f",big.mark=",",digits=0,width=11),
                   EAF=round(EAF,3),b=round(b,3),SE=round(SE,3),P=format(P,digits=3,scientific=TRUE),
                   Info=round(Info,3),
                   N=formatC(N,format="f",big.mark=",",digits=0,width=5)) %>%
            select(Protein,Chromosome,Position,rsid,EAF,b,SE,P,Info,Imputed)
 ukb_ppp <- read.delim(file.path(INF,"work","UKB-PPP.txt"))
    cojo <- merge(read.csv(file.path(INF,"sentinels","INF1.jma-rsid.cis.vs.trans.csv")),
                  gap_inf1[c("prot","target.short")],by="prot") %>%
            mutate(prot=target.short,log10P=round(-log10p(b/se),2),log10PJ=round(-log10p(bJ/bJ_se),2)) %>%
            rename(UniProt=uniprot,Protein=prot,Chromosome=Chr,Position=bp,rsid=SNP,SE=se,P=p,
                   bJ_SE=bJ_se,cistrans=cis.trans) %>%
            mutate(Position=formatC(Position,format="f",big.mark=",",digits=0,width=11),
                   b=round(b,3),SE=round(SE,3),bJ=round(bJ,3),bJ_SE=round(bJ_SE,3)) %>%
            select(UniProt,Protein,Chromosome,Position,rsid,b,SE,log10P,bJ,bJ_SE,log10PJ,cistrans)
   h2pve <- read.csv(file.path(INF,"ldak","h2-ldak-pve.csv"))
     vep <- merge(read.sheet("VEP", 1:27, 2:182),gap_inf1,by.x="Protein",by.y="prot") %>%
            mutate(Protein=target.short) %>% select(-target.short)
   eqtls <- read.sheet("eQTLs", 1:24, 2:24)
 eQTLGen <- read.table(file.path(INF,"eQTLGen","SCALLOP-INF-eQTLGen.txt"),header=TRUE) %>%
            left_join(gap_inf1) %>%
            mutate(prot=target.short,flag=if_else(P_eqtl<5e-8,"x","")) %>%
            rename(Protein=prot) %>% select(-target.short)
 eQTLGen_coloc <- read.table(file.path(INF,"eQTLGen","coloc.txt"),header=TRUE) %>%
                  rename(Gene=gene,H0=PP.H0.abf,H1=PP.H1.abf,H2=PP.H2.abf,H3=PP.H3.abf,H4=PP.H4.abf) %>%
                  left_join(gap_inf1) %>%
                  mutate(prot=target.short,
                         nSNP=nsnps,
                         H0=round(H0,2),
                         H1=round(H1,2),
                         H2=round(H2,2),
                         H3=round(H3,2),
                         H4=round(H4,2)) %>%
                  rename(Protein=prot,UniProt=uniprot) %>%
                  select(ID,UniProt,Protein,Gene,nSNP,H0,H1,H2,H3,H4)
 eQTLCatalogue <- read.delim(file.path(INF,"eQTLCatalogue","eQTLCatalogue-all.tsv"),header=TRUE) %>%
                  left_join(gap_inf1) %>%
                  mutate(prot=target.short,
                         nSNP=nsnps,
                         H0=round(H0,2),
                         H1=round(H1,2),
                         H2=round(H2,2),
                         H3=round(H3,2),
                         H4=round(H4,2)) %>%
                  rename(UniProt=uniprot,Protein=prot,SNPid=snpid,Study=unique_id) %>%
                  select(UniProt,Protein,rsid,Study,nSNP,H0,H1,H2,H3,H4)
reactome <- read.sheet("Reactome", 1:19, 2:589)
garfield <- read.table(file.path(INF,"garfield-data","output","INF1-cis","garfield.test.INF1.out"),header=TRUE) %>%
            rename(P=Pvalue,cellType=Celltype,b=Beta,LCL=CI95_lower,UCL=CI95_upper) %>%
            mutate(PThresh=format(PThresh,digits=3,scientific=TRUE),P=round(P,3),OR=round(OR,3),b=round(b,3),SE=round(SE,3)) %>%
            select(ID,PThresh,P,Annotation,cellType,Tissue,Type,Category,OR,b,SE)
   magma <- read.delim(file.path(INF,"work","All.dat"))
  fusion <- read.sheet("FUSION", 1:26, 2:117)
     smr <- merge(read.sheet("SMR", 1:27, 2:83),gap_inf1,by="prot") %>%
            mutate(prot=target.short) %>% rename(Protein=prot) %>% select(-target.short)
            d <- read.sheet("GSMR", 1:12, 2:55)
            na1 <- with(d,is.na(Exposure1))
            na2 <- with(d,is.na(Exposure2))
            d[na1,"Exposure1"] <- d[na1,"Exposure2"]
            d[na2,"Exposure2"] <- d[na2,"Exposure1"]
    gsmr <- merge(d, gap_inf1[c("prot","target.short")],by.x="Exposure1",by.y="prot") %>%
            mutate(Exposure1=target.short,Exposure2=target.short) %>% rename(Protein1=Exposure1,Protein2=Exposure2) %>%
            select(-target.short)
    gsmr_efo <- read.delim(file.path(INF,"mr","gsmr","gsmr-efo-reduce.txt")) %>%
                mutate(protein=gsub("MCP-1","CCL2",protein),
                       protein=gsub("MCP-2","CCL8",protein),
                       protein=gsub("MCP-3","CCL7",protein),
                       protein=gsub("MCP-4","CCL13",protein),
                       r2=NA) %>%
                rename(Protein=protein,ID=id,nSNP=nsnp,SE=se,
                       FDR=fdr,nCase=Ncase,nControl=Ncontrol,nTotal=Ntotal,pQTL=pqtl,P=p,QTL=qtl,P_QTL=p_qtl) %>%
                mutate(FDR=format(FDR,digits=3,scientific=TRUE),
                       nCase=formatC(nCase,format="f",big.mark=",",digits=0,width=11),
                       nControl=formatC(nControl,format="f",big.mark=",",digits=0,width=11),
                       nTotal=formatC(nTotal,format="f",big.mark=",",digits=0,width=11),
                       bxy=round(bxy,3),SE=round(SE,3),P=format(P,digits=3,scientific=TRUE),
                       P_QTL=format(P_QTL,digits=3,scientific=TRUE))
     for (i in 1:nrow(gsmr_efo))
     {
         z <- gsmr_efo[i,]
         r <- ieugwasr::ld_matrix_local(z[c("pQTL","QTL")],with_alleles=TRUE,
                                        bfile=file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",z$chr)),
                                        plink_bin="/rds/user/jhz22/hpc-work/bin/plink")
         r2 <- ifelse(nrow(r)==2,r[1,2]^2,r^2)
         gsmr_efo[i,"r2"] <- round(r2,3)
     }
     gsmr_efo <- select(gsmr_efo,-chr)
     crp <- read.sheet("CRP", 1:15, 2:30)
     gdb <- read.sheet("geneDrugbank", 1:7, 2:72)
     at1 <- readWorkbook(xlsxFile=url,sheet="Annotrans1"); #names(at1) <- replace(names(at1),grepl("^[X]",names(at1)),"")
     at2 <- readWorkbook(xlsxFile=url,sheet="Annotrans2"); #names(at2) <- replace(names(at2),grepl("^[X]",names(at2)),"")
     at3 <- readWorkbook(xlsxFile=url,sheet="Annotrans3"); #names(at3) <- replace(names(at3),grepl("^[X]",names(at3)),"")

  great3 <- read.delim(file.path(INF,"GREAT","IL12B-KITLG-TNFSF10.tsv")) %>%
            mutate(fdr=p.adjust(BinomP,method="fdr")) %>% arrange(fdr)
   great <- read.delim(file.path(INF,"GREAT","cistrans.tsv")) %>%
            mutate(fdr=p.adjust(BinomP,method="fdr")) %>% arrange(fdr)

read_table <- function(f, exprs="pval <= 0.05/nrow(t)")
{
  addflag <- function(exprs)
  {
    t <- within(read.delim(f), {flag=""})
    x <- with(t,eval(str2expression(exprs)))
    x <- with(t,eval(str2expression(e)))
    t[x, "flag"] <- "x"
  }
  t <- within(read.delim(f), {fdr <- p.adjust(pval,method="fdr")})
# t <- addflag(exprs)
  t
}

mr_immun <- merge(read_table(file.path(INF,"mr","pQTLs","pQTL-efo.txt")),gap_inf1,by.x="exposure",by.y="prot") %>%
            mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)
mr_misc <- merge(read_table(file.path(INF,"mr","pQTLs","pQTL-ieu-FEV1.txt")),gap_inf1,by.x="exposure",by.y="prot") %>%
           mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)

mr_read <- function(type)
{
  merge(read_table(file.path(INF,"mr",paste(type,"efo-result.txt",sep="-"))),
        gap_inf1, by.x="exposure",by.y="prot") %>%
  mutate(exposure=target.short) %>% rename(Protein=exposure) %>% select(-target.short) %>% arrange(fdr)
}
cis_mr <- mr_read("cis")
trans_mr <- mr_read("trans")
pan_mr <- mr_read("pan")
mr <- rbind(cis_mr,trans_mr,pan_mr)

protein_correlation <- read.csv(file.path(INF,"coffeeprot","table_complex.csv")) %>%
                       arrange(desc(cor)) %>%
                       mutate(varID1=toupper(varID1),varID2=toupper(varID2),VarVar=toupper(VarVar))
protein_dgi <- read.csv(file.path(INF,"coffeeprot","protein_annotated.csv")) %>%
               select(varID,ID,HPA_IF_protein_location,CP_loc,inDGIdb) %>%
               mutate(gene=toupper(ID)) %>%
               left_join(gap.datasets::inf1[c("gene","target.short")]) %>%
               rename(Protein=target.short) %>%
               select(-c(varID,ID,gene))
pqtl_annotation <- read.csv(file.path(INF,"coffeeprot","table_qtl_processed.csv")) %>%
                   rename(gene=gene_symbol) %>%
                   left_join(gap.datasets::inf1[c("gene","target.short")]) %>%
                   rename(Protein=target.short) %>%
                   select(-pvalue)
protein_dgi <- protein_dgi %>% select(Protein,names(protein_dgi))
pqtl_annotation <- pqtl_annotation %>% select(Protein,names(pqtl_annotation))
pqtl_impact <- read.csv(file.path(INF,"rsid","variant_effect_impact.csv"))

pav <- merge(within(pqtls,{prot_rsid=paste0(prot,"-",rsid)}),
             within(vep,{prot_rsid=paste0(Protein,"-",vep[["#Uploaded_variation"]])}),by="prot_rsid")
data.frame(table(subset(pav,cis.trans=="cis")$Consequence))

print(head(interval))
print(head(os))
print(head(cvd1))
print(head(fenland))
print(head(decode))
print(head(aric))
names(interval)[5] <- names(os)[5] <- names(cvd1)[5] <- names(fenland)[5] <- names(decode)[5] <- names(aric)[5] <- names(ages)[5] <- "cis/trans"
knownpqtls_dup <- bind_rows(interval,os,cvd1,fenland,decode,aric,ages,ngs) %>%
                  mutate(Protein=gsub("MCP-1","CCL2",Protein)) %>%
                  mutate(Protein=gsub("MCP-2","CCL8",Protein)) %>%
                  mutate(Protein=gsub("MCP-3","CCL7",Protein)) %>%
                  mutate(Protein=gsub("MCP-4","CCL13",Protein))
knownpqtls <- distinct(knownpqtls_dup[c("Sentinels","SNPid","UniProt","Protein")]) %>% arrange(Protein,SNPid)
pqtlstudies <- unique(knownpqtls_dup[c("Source","PMID")]) %>% arrange(PMID)
rownames(pqtlstudies) <- seq(nrow(pqtlstudies))

options("openxlsx.borderColour"="#4F80BD")
hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
url <- "https://jhz22.user.srcf.net/pqtl-immune_infection_edited.xlsx"
url <- file.path(INF,"work","pqtl-immune_infection_edited.xlsx")
credibleset <- read.table(file.path(INF,"work","INF1.merge-rsid.cs"),col.names=c("prot","MarkerName","CredibleSet"),sep="\t")
credibleppa <- read.table(file.path(INF,"work","INF1.merge-rsid.ppa"),col.names=c("prot","MarkerName","PPA"),sep="\t")
credibleset_unprune <- read.table(file.path(INF,"cs","unprune","INF1.merge-rsid.cs"),col.names=c("prot","MarkerName","CredibleSet_unpruned"),sep="\t")
credibleppa_unprune <- read.table(file.path(INF,"cs","unprune","INF1.merge-rsid.ppa"),col.names=c("prot","MarkerName","PPA_unpruned"),sep="\t")
pqtls <- merge(pqtls,credibleset,by.x=c("prot","rsid"),by.y=c("prot","MarkerName")) %>%
         merge(credibleppa,by.x=c("prot","rsid"),by.y=c("prot","MarkerName")) %>%
         merge(credibleset_unprune,by.x=c("prot","rsid"),by.y=c("prot","MarkerName")) %>%
         merge(credibleppa_unprune,by.x=c("prot","rsid"),by.y=c("prot","MarkerName")) %>%
         rename(Protein=prot,SNPid=MarkerName,cistrans=cis.trans,EAF=Freq1,b=Effect,SE=StdErr,logP=log.P.) %>%
         mutate(prots=Protein,
                UniProt=uniprot,
                Protein=target.short,
                Position=formatC(Position,format="f",big.mark=",",digits=0,width=11),
                EAF=round(EAF,digits=2),
                b=format(b,digits=2,justify="right"),
                SE=format(SE,digits=2,justify="right"),
                logP=format(-logP,digits=2,justify="right"),
                HetChiSq=formatC(HetChiSq,format="f",big.mark=",",digits=1,width=6),
                logHetP=format(-logHetP,digits=3,scientific=TRUE,justify="right"),
                N=formatC(N,format="f",big.mark=",",digits=0,width=5),
                lengthCS=unlist(lapply(sapply(CredibleSet,function(x) strsplit(x," ")),length)),
                lengthCS_unpruned=unlist(lapply(sapply(CredibleSet_unpruned,function(x) strsplit(x," ")),length))) %>%
         select(UniProt,Protein,Chromosome,Position,Start,End,cistrans,rsid,SNPid,Allele1,Allele2,EAF,
                b,SE,logP,Direction,HetISq,HetChiSq,HetDf,logHetP,N,
                CredibleSet,PPA,lengthCS,CredibleSet_unpruned,PPA_unpruned,lengthCS_unpruned,prots,uniprot,-target.short,-Start,-End)
metal <- read.delim(file.path(INF,"work","INF1.METAL"))
pqtldisease <- subset(read.sheet("short",1:51,1:220),Keep==1) %>%
               left_join(unique(pqtls[c("Protein","prots")]),by="prots") %>%
               mutate(prots=if_else(grepl(";",prots),prots,str_replace(prots,prots,Protein))) %>%
               mutate(prots=if_else(grepl("CXCL9;IL.12B",prots),str_replace_all(prots,c("IL.12B"="IL-12B")),prots)) %>%
               mutate(prots=if_else(grepl("MMP.10;CST5",prots),str_replace_all(prots,c("MMP.10"="MMP-10")),prots)) %>%
               rename(Proteins=prots) %>%
               left_join(distinct(metal[c("MarkerName","rsid")])) %>%
               rename(Trait=trait,EFO=efo,Study=study,PMID=pmid,Dataset=dataset) %>%
               select(rsid,Proteins,Allele1,Allele2,Effects,SEs,cistrans,Trait,EFO,Study,PMID,Dataset)
coloc <- merge(read.delim(file.path(INF,"coloc","GTEx-all.tsv")),gap_inf1,by="prot") %>%
         mutate(prot=target.short,
                nSNP=nsnps,
                H0=round(H0,2),
                H1=round(H1,2),
                H2=round(H2,2),
                H3=round(H3,2),
                H4=round(H4,2)) %>%
         rename(UniProt=uniprot,Protein=prot,SNPid=snpid,Tissue=qtl_id) %>%
         select(UniProt,Protein,rsid,Tissue,nSNP,H0,H1,H2,H3,H4)
cs95 <- read.delim(file.path(INF,"coloc-jma","cis-eQTL_table.tsv"))
cs95 <- data.frame(rsidProt=str_replace(rownames(cs95),"[.]","-"),cs95)
HOME <- Sys.getenv("HOME")
load(file.path(HOME,"software-notes","docs","files","pi_database.rda"))
drug <- subset(pi_drug,target%in%with(gap.datasets::inf1,gene)) %>% left_join(pi_trait) %>%
        transmute(Trait=trait,Target=target,Drug=drug,Mechanism=mechanism_of_action,actionType=action_type,Name=name,Source=source,maxPhase=max_phase)
efo <- read.delim(file.path(INF,"rsid","efo.txt"))
hgi_gsmr <- read.delim(file.path(INF,"mr","gsmr","hgi","5e-8","5e-8.tsv"))
hgi_pqtlmr <- read.delim(file.path(INF,"HGI","pqtlMR.txt"))
pqtls <- select(pqtls,-prots,-uniprot)
pqtls_o <- select(pqtls,-SNPid)

outsheets <- c("summary",
               "inf1","pqtls_o","cojo","aristotl",
                         "eQTLGen_coloc","coloc","eQTLCatalogue","pqtldisease","gsmr_efo","drug",
               "studies","garfield","pqtl_impact","vep","magma","hgi_gsmr","hgi_pqtlmr",
	       "knownpqtls","eQTLGen","reactome","great","efo","gdb",
               "interval","os","cvd1","fenland","decode","aric","ages","ukb_ppp","ngs","pqtlstudies",
               "great3","mr_immun","smr","cis_mr","mr_misc",
               "protein_correlation", "protein_dgi")
titles <- c("summary",
            "Inflammation-panel","pQTLs","Conditional-analysis","ARISTOTLE-study",
                      "eQTLGen-coloc","GTEx-coloc","eQTL-Catalogue-coloc","Disease-GWAS-overlap","GSMR-results","PI-drug",
            "Cohorts","GARFIELD-outputs","pQTL impact","VEP annotation","MAGMA outputs","HGI-GSMR r6","HGI-pQTLMR",
            "known pQTLs","eQTLGen","Reactome","GREAT","EFO","geneDrugbank",
            "INTERVAL study","Other studies","SCALLOP-CVD1","Fenland study","deCODE study","ARIC study","AGES study","UKB-PPP","UKB47k","previous pQTL studies",
            "IL12B-KITLG-TNFSF10","pQTL-immune-MR","SMR","cis-MR results","pQTL-misc-MR",
            "Protein correlation","DGI membership")
description=paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))
uppered <- c("PQTLs")
description[description%in%uppered] <- titles[description%in%uppered]
n0 <- 1
n1 <- 10
prefix <- c(paste0(toupper(substr(outsheets, 1, 1)), substr(outsheets, 2, nchar(outsheets)))[1:n0],
            rep("ST",n1),
            paste0(toupper(substr(titles, 1, 1)), substr(titles, 2, nchar(titles)))[(n0+n1+1):length(outsheets)]
          )
summary <- data.frame(Sheetnames=prefix,Description=description)[1:(n0+n1),]
summary2 <- data.frame(Sheetnames=prefix,Description=description)[-(1:(n0+n1)),]
xlsx <- file.path(INF,"NG","Supplementary-Tables.xlsx")
wb <- createWorkbook(xlsx)
xlsx2 <- file.path(INF,"NG","Additional-Tables.xlsx")
wb2 <- createWorkbook(xlsx2)
addWorksheet(wb2,"Summary",zoom=150)
writeData(wb2,"Summary","Summary",xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
          fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
writeDataTable(wb2, "Summary", summary2, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
for (i in 1:length(outsheets))
{
  if (i<=n0+n1)
  {
    sheetnames <- with(summary[i,], ifelse(i<=n0|i>n0+n1, Description, paste0(Sheetnames,"-",Description)))
    cat(sheetnames,"\n")
    addWorksheet(wb, sheetnames, zoom=150)
    writeData(wb, sheetnames, sheetnames, xy=c(1,1),
                  headerStyle=createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
    body <- get(outsheets[i])
    writeDataTable(wb, sheetnames, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
    freezePane(wb, sheetnames, firstCol=TRUE, firstActiveRow=3)
    width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
  # width_vec_header <- nchar(colnames(body))+2
    setColWidths(wb, sheetnames, cols = 1:ncol(body), widths = width_vec)
    writeData(wb, sheetnames, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
  } else {
    sheetnames <- with(summary2[i-n0-n1,], ifelse(i<=n0|i>n0+n1, Description, paste0(Sheetnames,"-",Description)))
    cat(sheetnames,"\n")
    addWorksheet(wb2, sheetnames, zoom=150)
    writeData(wb2, sheetnames, sheetnames, xy=c(1,1),
                  headerStyle=createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
    body <- get(outsheets[i])
    writeDataTable(wb2, sheetnames, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
    freezePane(wb2, sheetnames, firstCol=TRUE, firstActiveRow=3)
    width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
    setColWidths(wb2, sheetnames, cols = 1:ncol(body), widths = width_vec)
    writeData(wb2, sheetnames, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
 }
}
sheets_wb <- sheets(wb)
data.frame(sheets_wb)

bStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
hStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
#conditionalFormatting(wb, sheets_wb[grepl("IL12B-KITLG-TNFSF10",sheets_wb)], cols = 26, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GREAT",sheets_wb)], cols = 25, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GTEx coloc$",sheets_wb)], cols = 12, rows = 3:nrow(coloc), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("GARFIELD",sheets_wb)], cols = 3, rows = 3:nrow(garfield), rule = "<=1e-5", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("immune-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_immun), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("MR results",sheets_wb)], cols = 9, rows = 3:nrow(mr), rule = "==\"x\"", style = hStyle)
#conditionalFormatting(wb, sheets_wb[grepl("misc-MR",sheets_wb)], cols = 7, rows = 3:nrow(mr_misc), rule = "==\"x\"", style = hStyle)

saveWorkbook(wb, file=xlsx, overwrite=TRUE)
saveWorkbook(wb2, file=xlsx2, overwrite=TRUE)

# mr_immun <- read.sheet("pqtlMR-immune", 1:7, 2:67)
#  mr_misc <- read.sheet("pqtlMR-misc", 1:7, 2:39)
#      ivw <- read.sheet("IVW", 1:8, 2:19)
#     mrc2 <- read.sheet("MRC2", 1:7, 2:8)
#     mvmr <- read.sheet("MVMR", 1:9, 2:6)

annotate <- read.table(file.path(INF,"circos","annotate.txt"),header=TRUE) %>%
            mutate(gene=gsub("[.]"," ",gene))
novel_data <- subset(within(pqtls,{
                                    chrpos=paste0(Chromosome,":",Position)
                                    a1a2=paste0(toupper(Allele1),"/",toupper(Allele2))
                                    bse=paste0(b," (",SE,")")
                                    log10p=logP
                                  }),
                     !paste0(Protein,"-",rsid)%in%with(knownpqtls,paste0(Protein,"-",Sentinels)),
                     select=c(UniProt,Protein,SNPid,chrpos,rsid,a1a2,bse,log10p,cistrans,Chromosome,Position)) %>%
              left_join(annotate[c("uniprot","prot","p.gene","gene","MarkerName")],by=c('UniProt'='uniprot','SNPid'='MarkerName')) %>%
              rename(cis=cistrans,g.target=p.gene,g.pQTL=gene) %>%
              arrange(Chromosome,Position)
save(novel_data,file=file.path(INF,"work","novel_data.rda"))
prot_rsid <- with(novel_data,paste0(prot,"-",rsid))
prot_rsid_repl <- with(ukb_ppp,paste0(SCALLOP.prot,"-",SCALLOP.rsid))
left <- setdiff(prot_rsid,prot_rsid_repl)
if (FALSE)
{
  novel_data <- mutate(novel_data,prot_rsid=paste0(prot,"-",rsid)) %>%
                filter(prot_rsid %in% left) %>%
                select(-prot_rsid)
}
novelpqtls <- select(novel_data,Protein,chrpos,rsid,a1a2,bse,log10p,cis,g.target,g.pQTL,UniProt)
f <- file.path(INF,"NG","trans-pQTL_annotation.xlsx")
transpqtls <- read.xlsx(f,sheet=1,startRow = 6,colNames=FALSE,cols=c(1:3,14),skipEmptyRows=TRUE) %>%
              rename(Sentinel=X2,Encoding.gene=X3,causal.gene=X4) %>%
              mutate(causal.gene=if_else(causal.gene=="-","",causal.gene))
novelpqtls <- data.frame(no=1:nrow(novelpqtls),novelpqtls[,-c(1,10)]) %>%
              mutate(Sentinel=rsid) %>%
              left_join(transpqtls) %>%
              mutate(g.target=if_else(rsid=="rs3184504",Encoding.gene,g.target),
                     g.pQTL=if_else(is.na(causal.gene)|causal.gene=="",g.pQTL,paste0("*",causal.gene))) %>%
              select(-X1,-Sentinel,-Encoding.gene,-causal.gene)
write.xlsx(novelpqtls, file=file.path(INF,"NG","novelpqtls.xlsx"), overwrite=TRUE,
           colNames=TRUE,
           borders="surrounding", headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")

# known and novel pQTLs should be exclusive.
s <- pqtls %>% mutate(prot_rsid=paste0(Protein,"-",rsid))
s1 <- novel_data %>% mutate(prot_rsid=paste0(Protein,"-",rsid))
s2 <- knownpqtls %>% mutate(prot_rsid=paste0(Protein,"-",Sentinels))
c(s1$prot_rsid,s2$prot_rsid)[!c(s1$prot_rsid,s2$prot_rsid) %in% s$prot_rsid]

# ST-GARFIELD-outputs		This involves 1,005 combinations (42 Chromatin_States,  6 Formaldehyde-Assisted Isolation of Regulatory Elements (FAIREs), 41 Footprints. 7 Genic. 55 Histone_Modifications, 424 Hotspots, 424 Peaks, 6  transcription factor binding sites (TFBSs)) for four P value thresholds: 1e-2, 1e-5,1e-8, 1e-15, leading to enrichment P value with a Bonferroni cutoff of 0.05/4,020=1.24e-5. GARFIELD also has the garfield-Meff-Padj.R utility giving the effective number of annotations (457.04) and 1.1x10-4, leading to annotations reaching the association P=1.21x10-6 threshold with any protein.


#!/usr/bin/bash

function forestplot()
{
  cd ${INF}/rsid/10-4-2024
  module load ceuadmin/phenoscanner
  phenoscanner --snp=rs2228145 -c All -x EUR -r 0.8
# awk 'NR==1 || /CD40/' ${INF}/mr/gsmr/gsmr-efo-reduce-more.txt > gsmr-efo-reduce-CD40.txt
  R --no-save -q <<\ \ END
  options(width=200)
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
  metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
           filter(cis.trans=="cis") %>%
           select(prot,rsid,Chromosome) %>%
           rename(pqtl=rsid,chr=Chromosome)
  gsmr <- function (input="gsmr.txt",output="gsmr-efo.txt",top="gsmr-top.txt")
  {
    gsmr <- read.delim(input)
    efo <- read.delim(file.path(INF,"OpenGWAS","efo-update.txt")) %>%
           select(id,trait,ncase,ncontrol) %>%
           mutate(Ntotal=ncase+ncontrol) %>%
           filter(id %in% unique(gsmr$Outcome))
    gwas <- read.table(top,col.names=c("prot","id","qtl","a1_qtl","a2_qtl","freq","b_qtl","se_qtl","p_qtl","n_qtl")) %>%
            select(prot,id,qtl,b_qtl,se_qtl,p_qtl)
    gsmr_efo <- left_join(gsmr,pQTLdata::inf1[c("prot","target.short")], by=c("Exposure"="prot")) %>%
                left_join(filter(metal,prot %in% (gsmr$Exposure)), by=c("Exposure"="prot")) %>%
                left_join(efo,by=c("Outcome"="id")) %>%
                left_join(gwas,by=c("Exposure"="prot","Outcome"="id")) %>%
                rename(protein=Exposure,id=Outcome,Disease=trait,Ncase=ncase,Ncontrol=ncontrol) %>%
                mutate(protein=target.short,fdr=p.adjust(p,method="fdr")) %>%
                select(protein,Disease,id,nsnp,fdr,Ncase,Ncontrol,Ntotal,bxy,se,pqtl,p,qtl,b_qtl,se_qtl,p_qtl,chr,pqtl) %>%
                arrange(fdr)
    write.table(gsmr_efo,output,row.names=FALSE,quote=FALSE,sep="\t")
  }
  gsmr(input=file.path(INF,"mr","gsmr","gsmr-reduce.txt"),output="gsmr-efo-reduce-more.txt",top="gsmr-reduce-top.txt")
  CD40 <- read.delim("gsmr-efo-reduce-more.txt") %>%
          filter(protein=="CD40" &
                 id %in% c("ieu-a-833","ieu-b-18","ebi-a-GCST004131","ebi-a-GCST004132","ebi-a-GCST004133")) %>%
          select(Disease,bxy,se) %>%
          setNames(c("outcome","Effect","StdErr")) %>%
          mutate(outcome=gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))
  suppressMessages(library(gap))
  pdf("CD40.pdf",height=3,width=9)
  mr_forestplot(CD40,colgap.forest.left="0.05cm", fontsize=14,
                leftcols=c("studlab"), leftlabs=c("Outcome"),
                plotwidth="4inch", sm="OR",
                rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
                digits=2, digits.pval=2, scientific.pval=TRUE,
                common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
                addrow=TRUE, backtransf=TRUE, at=c(0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5), spacing=1.5, xlim=c(0.6,1.5))
  dev.off()
  # rs2228145, chr1:154426970, a1=A, a2=C
  traits <- "heumatoid|Coronary artery disease|abdominal aortic aneurysm|Abdominal aortic aneurysm|Atopic dermatitis"
  gwas <- read.delim("rs2228145_PhenoScanner_GWAS.tsv") %>%
          filter(grepl(traits,trait)) %>%
          filter(proxy==0 & r2==1 & !is.na(beta)) %>%
          select(-c(hg19_coordinates,ref_rsid,ref_hg19_coordinates,ref_a1,ref_a2,
                    ref_hg38_coordinates,hg38_coordinates,proxy,r2,dprime)) %>%
          rename(SNP=snp)
  swap <- gwas[["direction"]]=="-"
  gwas[swap,"a1"] <- "C"
  gwas[swap,"a2"] <- "A"
  write.table(gwas,file="ps.txt",sep="\t",quote=FALSE)
  select(gwas,study,efo,trait,beta,se,p,n)
  write.table(gwas[c(4,14,16),c("study","pmid","trait","efo","SNP","a1","a2","beta","se","p")],
              file="Ins.txt",sep="\t",row.names=FALSE,quote=FALSE)
  #########################################
  #MR pQTL analysis script -- Jie Zheng #
  #########################################
  #https://github.com/MRCIEU/epigraphdb-pqtl/blob/master/scripts/MR-pQTL-script.R

  ###install the Two sample MR package (just need do this once)
  ##source("https://bioconductor.org/biocLite.R")
  #install.packages("devtools")

  ##to update the R package (once there is a )
  #library(devtools)
  #install_github("MRCIEU/TwoSampleMR")

  #example of use the older version of the package
  #devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")

  ##call necessary libraries
  library(TwoSampleMR)
  #library(MRInstruments)
  #library("readxl") #plesae install the package
  #rm(list=ls(all=TRUE))

  ##setup your working folder
  #setwd("~/INF/rsid/jp-10-4-2024")

  ##read in the exposure data from a file
  #Ins<-data <- read_excel("./data/Instruments.xlsx",1)
  inf1 <- read.csv("INF1.merge.cis.vs.trans")
  subset(inf1,prot=="IL.6")
  jma <- read.delim("INF1.jma-rsid") %>%
         filter(prot=="IL.6") %>%
         mutate(Protein="IL-6",SNP="rs2228145",a1="A",a2="C")
  subset(jma,prot=="IL.6")
  Ins<-format_data(jma, type = "exposure", header = TRUE,
                   phenotype_col = "Protein", snp_col = "SNP", beta_col = "b",
                   se_col = "se", eaf_col = "freq", effect_allele_col = "a1",
                   other_allele_col = "a2", pval_col = "p")
  ##read in the outcome data
  #ao<-available_outcomes(access_token=NULL)

  #ids<-as.character(unlist(read.table("./data/outcome.id.txt",header=F)))

  #outcome_dat<- extract_outcome_data(
  #snps = Ins$SNP,
  #outcomes = ids)

  data <- read.delim("Ins.txt")
  outcome_dat <- format_data(data,
                             type = "outcome", header = TRUE,
                             phenotype_col = "trait", snp_col = "SNP", beta_col = "beta",
                             se_col = "se", eaf_col = "eaf", effect_allele_col = "a1",
                             other_allele_col = "a2", pval_col = "p")

  ##harmonise the exposure and outcome data
  dat <- NULL
  dat <- harmonise_data(exposure_dat = Ins, outcome_dat = outcome_dat)

  ##run the MR and sensitivity analyses
  mr_results <- NULL
  #mr_hetero <- NULL
  #mr_pleio <- NULL
  mr_single <- NULL
  #try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))  # main MR analysis
  try(mr_results <- mr(dat, method_list=c("mr_wald_ratio")))  # main MR analysis
  #mr_hetero <- mr_heterogeneity(dat) # heterogeneity test across instruments
  #mr_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test
  try(mr_single <- mr_singlesnp(dat)) #single SNP MR using Wald ratio

  ##save the MR results
  exposure <- "IL6-rs2228145"
  result_file0 <- paste0("./",exposure,".harmonise.txt")
  result_file <- paste0("./",exposure,".mr.txt")
  #result_file2 <- paste0("./",exposure,".mr_hetero.txt")
  #result_file3 <- paste0("./",exposure,".mr_pleio.txt")
  result_file4 <- paste0("./",exposure,".mr_single.txt")
  if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
  #if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
  #if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}
  IL6 <- mr_single %>%
         filter(grepl("^rs222",SNP)) %>%
         select(outcome,b,se) %>%
         setNames(c("outcome","Effect","StdErr")) %>%
         mutate(outcome=gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))
  library(gap)
  pdf("IL6.pdf",height=3,width=10)
  mr_forestplot(IL6,colgap.forest.left="0.05cm", fontsize=14,
                leftcols=c("studlab"), leftlabs=c("Outcome"),
                plotwidth="3.5inch", sm="OR",
                rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
                digits=2, digits.pval=2, scientific.pval=TRUE,
                common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
                addrow=TRUE, backtransf=TRUE, at=c(1:3)*0.5, spacing=1.5, xlim=c(0.5,1.5))
  dev.off()
  #as benchmark
  ins <- select(jma,prot,SNP,a1,a2,freq,b,se,p) %>%
         setNames(c("Phenotype","SNP","effect_allele","other_allele","eaf","beta","se","pval"))
  write.csv(ins,file="~/Ins.csv",quote=FALSE,row.names=FALSE)
  write.csv(rename(data,Phenotype=trait,effect_allele=a2,other_allele=a1,pval=p),
            file="~/Out.csv",quote=FALSE,row.names=FALSE)
  END
}

# 4,14,16
#  select(gwas,study,efo,trait,beta,se,n,direction)

#!/usr/bin/bash

Rscript -e '
  options(width=2000)
  require(dplyr)
  METAL <- read.delim("~/INF/work/INF1.METAL") %>%
           mutate(prot_rsid=paste0(uniprot,"-",rsid),start=if_else(Position-1e6<0,0,Position-1e6),end=Position+1e6)
  require(openxlsx)
  url <- "~/rds/results/public/proteomics/UKB-PPP/sun22.xlsx"
  ST10 <- read.xlsx(url,"ST10",startRow=3) %>%
          mutate(prot_rsid=paste0(Target.UniProt,"-",rsID))
  sentinels <- left_join(METAL,ST10,by="prot_rsid") %>%
               select(prot_rsid,cis.trans,rsID) %>%
               filter(!is.na(rsID))
  inf1 <- c(with(pQTLtools::inf1,uniprot),with(METAL,uniprot)) %>%
          unique()
  overlap <- filter(ST10,Target.UniProt %in% inf1)
  dim(overlap)
  UKB_PPP <- mutate(overlap,chrpos=strsplit(overlap[["Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)"]],":"),
                    chr=lapply(chrpos,"[[",1),
                    pos=lapply(chrpos,"[[",2),
                    chrpos=paste(chr,pos,sep=":"))
  suppressMessages(require(GenomicRanges))
  query <- with(UKB_PPP,GRanges(seqnames=as.integer(chr),IRanges(start=as.integer(pos),width=1),uniprot=Target.UniProt,rsid=rsID,prot=Assay.Target))
  subject <- with(METAL,GRanges(seqnames=Chromosome,IRanges(start=Position-1e6,width=2e6),uniprot=uniprot,rsid=rsid,pos=Position,prot=prot))
  fo <- findOverlaps(query,subject) %>% data.frame()
  co <- countOverlaps(query,subject) %>% data.frame()
  eq <- subject[fo$subjectHits,]$uniprot==query[fo$queryHits,]$uniprot
  ov <- data.frame(fo[eq,])
  ov1 <- query[ov$queryHits,] %>% data.frame() %>% select(-strand,-width) %>% mutate(rsid=if_else(rsid=="-",paste0("chr",seqnames,":",start),rsid))
  ov2 <- subject[ov$subjectHits,] %>% data.frame() %>% select(-strand,-width)
  b <- bind_cols(data.frame(ov1),data.frame(ov2))
  names(b) <- c(paste("UKB",names(ov1),sep="."),paste("SCALLOP",names(ov2),sep="."))
  token <- Sys.getenv("LDLINK_TOKEN")
  variant_list <- c(b[["UKB.rsid"]],b[["SCALLOP.rsid"]])
  r <- ieugwasr::ld_matrix(variant_list,pop="EUR",with_alleles=FALSE)
  excluded <- setdiff(variant_list,colnames(r))
  UKB.keep <- intersect(b$UKB.rsid,colnames(r))
  SCALLOP.keep <- intersect(b$SCALLOP.rsid,colnames(r))
  l <- matrix(NA,length(b$UKB.rsid),length(b$SCALLOP.rsid),dimnames=list(b$UKB.rsid, b$SCALLOP.rsid))
  l[UKB.keep,SCALLOP.keep] <- r[UKB.keep,SCALLOP.keep]
  r2 <-  sapply(1:nrow(b),function(x) with(b[x,],ifelse(UKB.rsid==SCALLOP.rsid,1,l[UKB.rsid,SCALLOP.rsid]^2)))
  replication <- mutate(b,r2=r2) %>%
                 filter(r2>=0.8)
  write.table(replication,file="~/INF/work/UKB-PPP.txt",row.names=FALSE,quote=FALSE,sep="\t")
  load("~/INF/work/novel_data.rda")
  prot_rsid <- with(novel_data,paste0(prot,"-",rsid))
  prot_rsid_repl <- with(replication,paste0(SCALLOP.prot,"-",SCALLOP.rsid))
  left <- setdiff(prot_rsid,prot_rsid_repl)

# variants need to be on the same chromosome
# r2 <- LDlinkR::LDmatrix(variant_list,pop="CEU",token=token)
# single pair
#  r <- sapply(1:nrow(b),function(x) with(b[x,],ifelse(UKB.rsid==SCALLOP.rsid,1,ieugwasr::ld_matrix(variants=c(UKB.rsid,SCALLOP.rsid),pop="EUR"))))
# r2 <- sapply(1:nrow(b),function(x) with(b[x,],ifelse(UKB.rsid==SCALLOP.rsid,1,LDlinkR::LDpair(UKB.rsid,SCALLOP.rsid,pop="CEU",token=token)$r2)))
'

Rscript -e '
  library(LDlinkR)
  LDexpress(snps = c("rs345", "rs456"),
                      pop = c("YRI", "CEU"),
                      tissue = c("ADI_SUB", "ADI_VIS_OME"),
                      r2d = "r2",
                      r2d_threshold = "0.1",
                      p_threshold = "0.1",
                      win_size = "500000",
                      genome_build = "grch37",
                      token = Sys.getenv("LDLINK_TOKEN")

'

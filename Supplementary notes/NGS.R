#
rds <- function(src="rds")
{
  if (src=="rds")
  {
    f <- "~/rds/results/public/proteomics/deCODE/doc/eldjarn22.xlsx"
    S1 <- read.xlsx(f,1)
    S3 <- read.xlsx(f,3)
  } else {
    f <- 'https://www.biorxiv.org/content/biorxiv/early/2022/02/20/2022.02.18.481034/DC3/embed/media-3.xlsx'
    S1 <- read.xlsx(f)
    f <- 'https://www.biorxiv.org/content/biorxiv/early/2022/02/20/2022.02.18.481034/DC12/embed/media-12.xlsx'
    S3 <- read.xlsx(f)
  }
  save(S1,S3,file=file.path(INF,"ukb","tables.rda"))
}

suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))

INF <- Sys.getenv("INF")
rds()
load(file.path(INF,"ukb","tables.rda"))
metal <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
         mutate(uniprot_rsid=paste0(uniprot,"-",rsid))
ss <- merge(S1,S3) %>%
      mutate(uniprot_rsid=paste0(uniprot,"-",rsName)) %>%
      filter(uniprot_rsid%in%metal$uniprot_rsid) %>%
      left_join(metal) %>%
      left_join(select(gap.datasets::inf1,prot,target,target.short)) %>%
      select(target.short,rsid,MarkerName,uniprot,cis.trans,Log10.pval.gc.cor.unadj) %>%
      transmute(Protein=target.short,Sentinels=rsid,UniProt=uniprot,SNPid=MarkerName,"cis/trans"=cis.trans,
                Proxies="as sentinel",r2=1,p=10^-Log10.pval.gc.cor.unadj,Source="Eldjarn, et al. (2022)",PMID="",Comment="")
options(width=200)
write.table(ss,file=file.path(INF,"ukb","NGS.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
cojo <- merge(read.csv(file.path(INF,"sentinels","INF1.jma-rsid.cis.vs.trans.csv")),gap.datasets::inf1[c("prot","target.short")],by="prot") %>%
        mutate(prot=target.short,uniprot_rsid=paste0(uniprot,"-",SNP)) %>% rename(Protein=prot) %>% select(-target.short)
region <- merge(S1,S3) %>%
          filter(uniprot %in% cojo[["uniprot"]])
key <- Sys.getenv("LDLINK_TOKEN")

signals <- function(src="metal")
{
  if (src=="metal") sentinels <- with(metal,uniprot_rsid)
  else sentinels <- with(subset(cojo,!uniprot_rsid %in% with(metal,uniprot_rsid)),uniprot_rsid)
  blocks <- r <- list()
  for(i in 1:length(sentinels))
  {
    if(src=="metal")
    {
      x <- subset(metal,uniprot_rsid==sentinels[i])
      snp <- with(x,rsid)
      prot <- with(x,prot)
      uniprot <- with(x,uniprot)
      chr <- with(x,paste0("chr",Chromosome))
      pos <- with(x,Position)
    } else {
      x <- subset(cojo,uniprot_rsid==sentinels[i])
      snp <- with(x,SNP)
      prot <- with(x,p.prot)
      uniprot <- with(x,uniprot)
      chr <- with(x,paste0("chr",Chr))
      pos <- with(x,bp)
    }
    blocks[[i]] <- subset(region,UniProt==uniprot & Chrom==chr & Pos>=pos-1e6 & Pos<pos+1e6)
    if(nrow(blocks[[i]])==0) next
    snps <- blocks[[i]][["rsName"]]
    sentinel_and_snps <- c(snp,snps[grepl("^rs",snps)])
    len <- length(sentinel_and_snps)
    cat("No", i,"prot-uniprot-pQTL", paste0(prot,"-",sentinels[i]),"chr =", chr, "pos =", pos, "total SNPs =", len, "\n")
    print(subset(select(blocks[[i]],Chrom,Pos,uniprot,GeneName,panel,rsName,CisOrTrans,Log10.pval.gc.cor.unadj),10^-Log10.pval.gc.cor.unadj<=5e-8))
    if(len >=2 & len <1000)
    {
       r[[i]] <- LDlinkR::LDmatrix(snps=sentinel_and_snps,pop="EUR",r2d="r2",token=key)
       r2 <- subset(r[[i]],RS_number==snp)
       sel <- !is.na(r2) & r2>=0.8
       print(names(r2)[sel])
       print(r2[sel])
    }
  }
}

signals("metal")
signals("cojo")

#!/usr/bin/bash

if [ ! -d ${INF}/GTEx ]; then mkdir ${INF}/GTEx; fi

# SNP annotation
phenoscanner -s T -c None -x EUR -p 5e-8 -r 0.8 -i ${INF}/sentinels/INF1.jma-rsid.cis -o ${INF}/GTEx/INF1.jma

# eQTL lookup
phenoscanner -s T -c eQTL -x EUR -p 5e-8 -r 0.8 -i ${INF}/sentinels/INF1.jma-rsid.cis -o ${INF}/GTEx/INF1.jma-rsid

R --no-save -q <<END
  library(dplyr)
  INF <- Sys.getenv("INF")
# annotation results
  jma <- read.table(file.path(INF,"sentinels","INF1.jma-rsid.cis.vs.trans"),as.is=TRUE,header=TRUE)
  INF1 <- within(left_join(subset(jma,cis.trans=="cis"),subset(gap.datasets::inf1,select=-c(start,end))),{
                 hg19_coordinates <- paste0("chr",Chr,":",bp)
                 HLA <- as.numeric(Chr==6 & bp >= 25392021 & bp <= 33392022)
          }) %>% rename(INF1_rsid=SNP, gene_gwas=gene, uniprot_gwas=uniprot)
  snps <- read.delim("GTEx/INF1.jma_PhenoScanner_SNP_Info.tsv") %>%
          mutate(snpid=gap::chr_pos_a1_a2(chr,pos_hg19,a1,a2)) %>%
          rename(gene=hgnc)
  INF1_aggr <- within(merge(INF1,snps,by.x="INF1_rsid",by.y="snp"), {gene_snpid <- paste0(gene,"-",snpid)})
# lookup results
  results <- read.delim(file.path(INF,"GTEx","INF1.jma-rsid_PhenoScanner_eQTL.tsv"))
  snps <- read.delim(file.path(INF,"GTEx","INF1.jma-rsid_PhenoScanner_SNP_Info.tsv")) %>%
          mutate(snpid=gap::chr_pos_a1_a2(chr,pos_hg19,a1,a2))
  ps <- right_join(snps[c("snp","rsid","hg19_coordinates","a1","a2","consequence","hgnc","proxy","r2","snpid")],
                   subset(results,study=="GTEx",select=-c(ref_rsid, ref_hg19_coordinates, ref_hg38_coordinates, dprime,
                                                          ref_a1, ref_a2, hg38_coordinates, efo, trait, probe, exp_ensembl,
                                                          ancestry, year, direction,n,n_studies,unit,dataset)))
  save(ps,file=file.path(INF,"GTEx","cis-pQTL.eQTL"))
END

# phenoscanner -t T -c eQTL -x EUR -p 5e-8 -r 0.8 -i work/INF1.gene -o INF1
# additionally from pqtlxQTL.R

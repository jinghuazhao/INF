#!/usr/bin/bash

R --no-save <<END
  library(gap)
  INF <- Sys.getenv("INF")
  INF1_prot <- read.table(file.path(INF,"work","INF1.merge.prot"),col.names=c("prot","uniprot"))
  prot <- subset(inf1[c("prot","target.short")], prot %in% with(INF1_prot,prot))
  write.table(prot,file=file.path(INF,"work","INF1.merge.target"),quote=FALSE,col.names=FALSE,row.names=FALSE)
END

# rsID,trait_name,snp_bp_abs,snp_chr,empirical_pvalue,trait_category
(
  echo rsID,trait_name,snp_chr,snp_bp_abs,LOD,trait_category
  sed '1d' ${INF}/work/INF1.METAL | \
  awk '{print $2,$3,$4,$5,-$11,$21}' | \
  sort -k2,2 | \
  join -22 -o2.1,1.2,2.3,2.4,2.5,2.6 <(sort -k1,1 ${INF}/work/INF1.merge.target) - | \
  tr ' ' ','
) > ${INF}/coffeeprot/qtl.csv

# rsID,snp_bp_abs,snp_chr,gene_symbol,gene_start_abs,gene_end_abs,gene_chr,pvalue,proxy
(
  echo rsID,snp_bp_abs,snp_chr,gene_symbol,gene_start_abs,gene_end_abs,gene_chr,LOD,proxy
  awk -v FS="," 'NR>1{print $5,$4,$3,$10,$12,$13,$11,-$6,$14}' ${INF}/work/INF1.merge.cis.vs.trans | \
  sort -k1,1 | \
  join ${INF}/work/INF1.merge.rsid - | \
  cut -d' ' -f1 --complement --output-delimiter=","
) > ${INF}/coffeeprot/pqtl.csv

head -1 work/INF1.merge.cis.vs.trans | tr ',' '\n' | awk '{print "#" NR,$1}'
#1 uniprot
#2 prot
#3 Chr
#4 bp
#5 SNP
#6 log10p
#7 p.prot
#8 p.target
#9 p.target.short
#10 p.gene
#11 p.chr
#12 p.start
#13 p.end
#14 cis.trans
#15 cis
#16 cis.end
#17 cis.start

(
  echo "chr","start_bp_mm10","end_bp_mm10"
  sed '1d;s/chr//' ${INF}/tryggve/EURLD.bed | \
  cut -f1-3 --output-delimiter=","
) > ${INF}/coffeeprot/EUR-LD.csv

R --no-save <<END
  require(Biobase)
  x <- readRDS("/home/jhz22/rds/post_qc_data/interval/phenotype/olink_proteomics/post-qc/eset.inf1.flag.out.outlier.in.rds")
  fn <- featureNames(x)
  m <- exprs(x)
  INF <- Sys.getenv("INF")
  a <- gap::inf1[c("prot","target.short")]
  b <- data.frame(prot=fn,m)
  INF1_prot <- read.table(file.path(INF,"work","INF1.merge.prot"),col.names=c("prot","uniprot"))
  prot <- subset(merge(a,b,by="prot"),prot%in%with(INF1_prot,prot),select=-prot)
  write.table(prot,
              file=file.path(INF,"coffeeprot","interval.csv"),quote=FALSE,row.names=FALSE,sep=",")
END

for suffix in ct intra ve vi
do
  convert \( mhtplot-${suffix}.png -append edges-${suffix}.png -append proteins-${suffix}.png -append \) +append mep-${suffix}.png
done

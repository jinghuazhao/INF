#!/usr/bin/bash

# cis-pQTLs
R --no-save <<END
  options(width=200)
  library(pQTLtools)
  INF <- Sys.getenv("INF")
  cvt <- read.csv(paste0(INF,"/work/INF1.merge.cis.vs.trans"))
  cvt_cis <- merge(subset(cvt,cis),inf1[c("gene","ensembl_gene_id")],by.x="p.gene",by.y="gene")
  cis_dat <- within(cvt_cis,{seqnames <- paste0("chr",Chr); start <- as.integer(bp-1); end <- as.integer(bp)})
  ord_cis <- with(cis_dat,order(Chr,bp))
  cis_dat <- data.frame(cis_dat[ord_cis,c("seqnames","start","end","SNP","uniprot","prot","p.gene","ensembl_gene_id","p.start","p.end")])
  cis_gr <- with(cis_dat,GenomicRanges::GRanges(seqnames=seqnames,IRanges::IRanges(start,end)))
  hpc_work <- Sys.getenv("HPC_WORK")
  path = file.path(hpc_work, "bin", "hg19ToHg38.over.chain")
  library(rtracklayer)
  chain <- import.chain(path)
  chain
  seqlevelsStyle(cis_gr) <- "UCSC"
  cis38 <- liftOver(cis_gr, chain)
  class(cis38)
  cis_dat38 <- data.frame(cis38)[c("seqnames","start","end")]
  names(cis_dat38) <- c("chr38","start38","end38")
  write.table(data.frame(cis_dat,cis_dat38),file=file.path(INF,"coloc","cis-pQTL.dat"),quote=FALSE,row.names=FALSE,sep="\t")
END

# cis-eQTLs from GTEx v8
export GTEx_v8=~/rds/public_databases/GTEx/GTEx_Analysis_v8_eQTL
export ext=.v8.signif_variant_gene_pairs.txt.gz
export col_gene=2
export col_variant=1
export GTEx_v8=~/rds/public_databases/GTEx/GTEx_Analysis_v8_eQTL_cis_associations
export ext=.v8.EUR.signif_pairs.txt.gz
export col_gene=1
export col_variant=2
export M=1e6
export nlines=60
(
  parallel -C' ' --env GTEx_v8 --env ext --env col --env M '
    read SNP hgnc ensGene pos chrpos < <(awk -v row={1} "NR==row{print \$4,\$7,\$8,\$13,\$11\"_\"\$13}" ${INF}/coloc/cis-pQTL.dat)
    zgrep ${ensGene} ${GTEx_v8}/{2}${ext} | \
    awk -v col_gene=${col_gene} -v col_variant=${col_variant} -v ensGene=${ensGene} -v M=${M} -v bp=${pos} -v OFS="\t" "
    {
       split(\$col_variant,a,\"_\");
       chr=a[1];pos=a[2];a1=a[3];a2=a[4];
       if(a[3]<a[4]) {a1=a[3];a2=a[4]} else {a1=a[4];a2=a[3]}
       \$col_variant=chr\":\"pos\"_\"a1\"_\"a2;
       snpid=chr\":\"bp\"_\"a1\"_\"a2
       if(index(\$col_gene,ensGene) && bp>=a[2]-M && bp<a[2]+M && length(a1)==1 && length(a2)==1) print \$col_gene,\$col_variant,\$7,chr,pos,a1,a2,snpid
    }" > ${INF}/coloc/${SNP}-${ensGene}-{2}.dat
    cat ${INF}/coloc/${SNP}-${ensGene}-{2}.dat | \
    sort -k3,3g | \
    awk -vSNP=${SNP} -vensGene=${ensGene} -vtissue={2} -vOFS="\t" "NR==1 {print SNP,ensGene,tissue,\$0}"
  ' ::: $(seq 2 ${nlines}) ::: $(ls ${GTEx_v8} | grep -v egenes | xargs -l basename -s ${ext})
) | \
find . -type f -empty -delete
sort -k1,1 -k2,2 > ${INF}/coloc/eQTL_GTEx.dat
awk '$5==$11' ${INF}/coloc/eQTL_GTEx.dat | cut -f1 | uniq
awk '$5==$11' ${INF}/coloc/eQTL_GTEx.dat | cut -f1 | uniq | grep -f - ${INF}/work/INF1.METAL | cut -f2,3
awk '$5==$11' ${INF}/coloc/eQTL_GTEx.dat

# Variants cis-eQTL regions, subject to check on r2
R --no-save <<END
  options(width=200)
  library(pQTLtools)
  eQTL_GTEx <- read.table("eQTL_GTEx.dat",col.names=c("MarkerName","ensGene","tissue","gene","sentinel","p","chr","pos","a1","a2","eQTL"))
  cis_dat <- within(eQTL_GTEx,{seqnames <- chr; start <- pos-1; end <- pos})
  ord_cis <- with(cis_dat,order(chr,pos))
  cis_dat <- cis_dat[ord_cis,c("seqnames","start","end","MarkerName","ensGene","tissue","gene","sentinel","eQTL","a1","a2")]
  cis_gr <- with(cis_dat,GenomicRanges::GRanges(seqnames=seqnames,IRanges::IRanges(start,end)))
  hpc_work <- Sys.getenv("HPC_WORK")
  path = file.path(hpc_work, "bin", "hg38ToHg19.over.chain")
  library(rtracklayer)
  chain <- import.chain(path)
  chain
  seqlevelsStyle(cis_gr) <- "UCSC"
  cis37 <- liftOver(cis_gr, chain)
  class(cis37)
  cis_dat37 <- data.frame(cis37)[c("seqnames","start","end")]
  names(cis_dat37) <- c("chr37","start37","end37")
  cisQTL <- within(cbind(cis_dat,cis_dat37),{snpid <- paste0(chr37,":",end37,"_",a1,"_",a2)})
  write.table(cisQTL,file="cis-eQTL.dat",quote=FALSE,row.names=FALSE,sep="\t")
END
for SNP in $(cut -f4 cis-eQTL.dat | sed '1d' | uniq)
do
  echo ${SNP}
  (
    echo ${SNP}
    awk -vSNP=${SNP} '$4==SNP' cis-eQTL.dat | cut -f15 | uniq
  ) > ${SNP}.snps
  plink --bfile INTERVAL/cardio/INTERVAL --extract ${SNP}.snps --make-bed --out ${SNP}
  plink --bfile ${SNP} --no-sex --no-pheno --r2 inter-chr --out ${SNP}
done
ls *.ld | xargs -I {} basename {} .ld | parallel -C' ' 'grep {} {}.ld'

# 95%CS
export cs95=~/rds/public_databases/GTEx/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz
export cs95tissue=${INF}/coloc/cs95tissue.dat
echo tissue> ${cs95tissue}
(
  sed '1d' ${INF}/coloc/cis-pQTL.dat | cut -f4,7,8,11,13 | awk '{print $1,$2,$3,$4"_"$5}' | \
  parallel -C' ' --env cs95 --env cs95tissue '
    zgrep {3} ${CS95} | zgrep {4} | \
    awk -v snpid={1} -v hgnc_symbol={2} -v ensGene={3} -v chrpos={4} -v OFS="\t" -v cs95tissue=${cs95tissue} "
    {
      printf snpid OFS hgnc_symbol OFS ensGene OFS chrpos OFS \$3
      split(\$6,a,\"|\")
      for(v in a) if(index(a[v],ensGene))
      {
        split(a[v],b,\"@\"); printf \" \" b[2];
        split(b[2],c,\"=\"); print c[1] >> cs95tissue
      }
      printf \"\n\"
    }"
  '
) | sort -k1,1 | join -t$'\t' <(cat ${INF}/work/INTERVAL.rsid | tr ' ' '\t') - > ${INF}/coloc/cis-eQTL-cs95.tsv

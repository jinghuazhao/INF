#!/usr/bin/bash

if [ ! -d ${INF}/coloc-jma ]; then mkdir ${INF}/coloc-jma; fi

# cis-pQTLs
R --no-save <<END
  options(width=200)
  library(pQTLtools)
  INF <- Sys.getenv("INF")
  cvt <- read.table(file.path(INF,"sentinels","INF1.jma-rsid.cis.vs.trans"),header=TRUE)
  cvt_cis <- merge(subset(cvt,cis),inf1[c("gene","ensembl_gene_id")],by.x="p.gene",by.y="gene")
  cis_dat <- within(cvt_cis,{seqnames <- paste0("chr",Chr); start <- as.integer(bp-1); end <- as.integer(bp)})
  ord_cis <- with(cis_dat,order(Chr,bp))
  cis_dat <- data.frame(cis_dat[ord_cis,c("seqnames","start","end","SNP","uniprot","prot","p.gene","ensembl_gene_id","p.start","p.end")])
  cis_gr <- with(cis_dat,GenomicRanges::GRanges(seqnames=seqnames,IRanges::IRanges(start,end)))
  f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
  library(rtracklayer)
  chain <- import.chain(f)
  chain
  seqlevelsStyle(cis_gr) <- "UCSC"
  cis38 <- liftOver(cis_gr, chain)
  class(cis38)
  cis_dat38 <- data.frame(cis38)[c("seqnames","start","end")]
  names(cis_dat38) <- c("chr38","start38","end38")
  write.table(data.frame(cis_dat,cis_dat38),file=file.path(INF,"coloc-jma","cis-pQTL.dat"),quote=FALSE,row.names=FALSE,sep="\t")
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
export nlines=100
parallel -C' ' --env GTEx_v8 --env ext --env col --env M '
    read SNP hgnc ensGene pos chrpos < <(awk -v row={1} "NR==row{print \$4,\$7,\$8,\$13,\$11\"_\"\$13}" ${INF}/coloc-jma/cis-pQTL.dat)
    zgrep ${ensGene} ${GTEx_v8}/{2}${ext} | \
    awk -v col_gene=${col_gene} -v col_variant=${col_variant} -v ensGene=${ensGene} -v M=${M} -v bp=${pos} -v OFS="\t" "
    {
       split(\$col_variant,a,\"_\");
       chr=a[1];pos=a[2];a1=a[3];a2=a[4];
       if(a[3]<a[4]) {a1=a[3];a2=a[4]} else {a1=a[4];a2=a[3]}
       \$col_variant=chr\":\"pos\"_\"a1\"_\"a2;
       snpid=chr\":\"bp\"_\"a1\"_\"a2
       if(index(\$col_gene,ensGene) && bp>=a[2]-M && bp<a[2]+M && length(a1)==1 && length(a2)==1) print \$col_gene,\$col_variant,\$7,chr,pos,a1,a2,snpid
    }" > ${INF}/coloc-jma/${SNP}-${ensGene}-{2}.dat
' ::: $(seq 2 ${nlines}) ::: $(ls ${GTEx_v8} | grep -v egenes | xargs -l basename -s ${ext})
(
  parallel -C' ' --env GTEx_v8 --env ext --env col --env M '
    read SNP hgnc ensGene pos chrpos < <(awk -v row={1} "NR==row{print \$4,\$7,\$8,\$13,\$11\"_\"\$13}" ${INF}/coloc-jma/cis-pQTL.dat)
    cat ${INF}/coloc-jma/${SNP}-${ensGene}-{2}.dat | \
    sort -k3,3g | \
    awk -vSNP=${SNP} -vensGene=${ensGene} -vtissue={2} -vOFS="\t" "NR==1 {print SNP,ensGene,tissue,\$0}"
  ' ::: $(seq 2 ${nlines}) ::: $(ls ${GTEx_v8} | grep -v egenes | xargs -l basename -s ${ext})
) | \
sort -k1,1 -k2,2 > ${INF}/coloc-jma/eQTL_GTEx.dat
awk '$5==$11' ${INF}/coloc-jma/eQTL_GTEx.dat | cut -f1 | uniq
awk '$5==$11' ${INF}/coloc-jma/eQTL_GTEx.dat | cut -f1 | uniq | grep -f - ${INF}/work/INF1.METAL | cut -f2,3
awk '$5==$11' ${INF}/coloc-jma/eQTL_GTEx.dat
find . -type f -empty -delete

# Variants cis-eQTL regions, subject to check on r2
R --no-save -q <<END
  options(width=200)
  library(pQTLtools)
  INF <- Sys.getenv("INF")
  eQTL_GTEx <- read.table(file.path(INF,"coloc-jma","eQTL_GTEx.dat"),
                          col.names=c("MarkerName","ensGene","tissue","gene","sentinel","p","chr","pos","a1","a2","eQTL"))
  cis_dat <- within(eQTL_GTEx,{seqnames <- chr; start <- pos-1; end <- pos})
  ord_cis <- with(cis_dat,order(chr,pos))
  cis_dat <- cis_dat[ord_cis,c("seqnames","start","end","MarkerName","ensGene","tissue","gene","sentinel","eQTL","a1","a2")]
  cis_gr <- with(cis_dat,GenomicRanges::GRanges(seqnames=seqnames,IRanges::IRanges(start,end)))
  library(rtracklayer)
  f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
  chain <- import.chain(f)
  chain
  seqlevelsStyle(cis_gr) <- "UCSC"
  cis37 <- liftOver(cis_gr, chain)
  class(cis37)
  cis_dat37 <- data.frame(cis37)[c("seqnames","start","end")]
  names(cis_dat37) <- c("chr37","start37","end37")
  cisQTL <- within(cbind(cis_dat,cis_dat37),{snpid <- paste0(chr37,":",end37,"_",a1,"_",a2)})
  write.table(cisQTL,file=file.path(INF,"coloc-jma","cis-eQTL.dat"),quote=FALSE,row.names=FALSE,sep="\t")
END
for SNP in $(cut -f4 ${INF}/coloc-jma/cis-eQTL.dat | sed '1d' | uniq)
do
  echo ${SNP}
  (
    echo ${SNP}
    awk -vSNP=${SNP} '$4==SNP' ${INF}/coloc-jma/cis-eQTL.dat | cut -f15 | uniq
  ) > ${INF}/coloc-jma/${SNP}.snps
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/coloc-jma/${SNP}.snps --make-bed --out ${INF}/coloc-jma/${SNP}
  plink --bfile ${INF}/coloc-jma/${SNP} --no-sex --no-pheno --r2 inter-chr --out ${INF}/coloc-jma/${SNP}
done
ls *.ld | xargs -I {} basename {} .ld | parallel -C' ' 'grep {} ${INF}/coloc-jma/{}.ld'

# 95%CS -- SNP already in rsid so no need to copy from INTERVAL.rsid
export DAPG=~/rds/public_databases/GTEx/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz
export cs95=${INF}/coloc-jma/cis-eQTL_cs95.tsv
export cs95tissue=${INF}/coloc-jma/cs95tissue.dat
echo tissue > ${cs95tissue}
(
  sed '1d' ${INF}/coloc-jma/cis-pQTL.dat | cut -f4,6,8,11,13 | awk '{print $1,$2,$3,$4"_"$5}' | \
  parallel -C' ' --env DAPG --env cs95tissue '
    zgrep {3} ${DAPG} | zgrep {4} | \
    awk -v snpid={1} -v hgnc_symbol={2} -v ensGene={3} -v chrpos={4} -v OFS="\t" -v cs95tissue=${cs95tissue} "
    {
      printf snpid OFS hgnc_symbol OFS ensGene OFS chrpos OFS \$3 OFS
      split(\$6,a,\"|\")
      for(v in a) if(index(a[v],ensGene))
      {
        split(a[v],b,\"@\"); printf \" \" b[2];
        split(b[2],c,\"=\"); print c[1] >> cs95tissue
      }
      printf \"\n\"
    }"
  '
) > ${cs95}
R --no-save <<END
  library(dplyr)
  eqtl_file <- Sys.getenv("cs95")
  eqtls <- read.table(eqtl_file,sep="\t", col.names=c("rsid","prot","ensGene","chrpos","GTExSNP","tissue_p")) %>%
           left_join(gap.datasets::inf1[c("prot","target.short")]) %>%
           mutate(rsidProt=paste0(rsid," (",target.short,")"),tissue_p=sub("^ ","",tissue_p)) %>%
           arrange(target.short)
  tissue_file <- Sys.getenv("cs95tissue")
  tissues <- with(read.table(tissue_file,header=TRUE),sort(unique(tissue)))
  eqtl_table <- matrix("",nrow(eqtls),length(tissues))
  rownames(eqtl_table) <- with(eqtls,rsidProt)
  colnames(eqtl_table) <- tissues
  for(row in with(eqtls,rsidProt))
  {
    all_pairs <- unlist(strsplit(subset(eqtls,rsidProt==row)[["tissue_p"]]," "))
    cat(row,length(all_pairs))
    for(tp in all_pairs)
    {
      z <- unlist(strsplit(tp,"="))
      tissue <- z[1]
      p <- z[2]
      cat(",",tissue)
      eqtl_table[row,tissue] <- gsub("^;","",paste(eqtl_table[row,tissue],p,sep=";"))
    }
    cat("\n")
  }
  library(stringr)
  colnames(eqtl_table) <- gsub("ba","BA",gsub("_"," ",str_to_sentence(colnames(eqtl_table))))
  INF <- Sys.getenv("INF")
  write.table(eqtl_table,file=file.path(INF,"coloc-jma","cis-eQTL_table.tsv"),sep="\t")
  tbl <- eqtl_table
  tbl[eqtl_table==""] <- 0
  tbl[eqtl_table!=""] <- 1
  storage.mode(tbl) <- "integer"
  library(pheatmap)
  pal <- colorRampPalette(c("white","red"))
  col <- pal(3)
  library(grid)
  png(file.path(INF,"coloc-jma","cis_eQTL.png"),res=300,width=15,height=10,units="in")
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(tbl, legend=FALSE, angle_col="315", color=col, width=12, height=40, cluster_rows=FALSE, cluster_cols=FALSE, fontsize=18)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Tissue", y=-0.07, gp=gpar(fontsize=15))
  grid.text("pQTL", x=-0.07, rot=90, gp=gpar(fontsize=15))
  dev.off()
END

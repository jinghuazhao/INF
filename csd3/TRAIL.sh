#!/usr/bin/bash

export prot=TRAIL
export uniprot=P50591

function turboman()
{
  (
    echo \#CHROM POS P
    gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
    cut -f1,2,12 | \
    sed '1d;s|\t| |g;s|-||g' | \
    sort -k1,1n -k2,2n
  ) | \
  gzip -f > ${INF}/work/${prot}.dat.gz

  (
    echo chromosome position
    grep ${prot} ${INF}/work/INF1.merge | \
    cut -f8,9 | \
    sort -k1,1n -k2,2n
  ) | \
  paste - <(echo SNP CFHR1\;CFHR3 TNFSF10 KNG1 MYRF SERPINA1 APOH MEP1B PLAUR | tr ' ' '\n') > ${INF}/work/${prot}.annotate

  export refgene_file_name=turboman_hg19_reference_data.rda
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    refgene_file_name <- Sys.getenv("refgene_file_name")
    load(file.path(INF,"cardio",refgene_file_name))
    library(dplyr)
    save(ld_block_breaks_pickrell_hg19_eur,refgene_gene_coordinates_h19,file=refgene_file_name,compress="xz")
  END

  R --slave --vanilla --args \
    input_data_path=${INF}/work/${prot}.dat.gz \
    output_data_rootname=${INF}/work/${prot} \
    custom_peak_annotation_file_path=${INF}/work/${prot}.annotate \
    reference_file_path=${refgene_file_name} \
    pvalue_sign=5e-10 \
    plot_title="${prot} (${uniprot})" < ${INF}/cardio/turboman.r
}

# turboman

# grep TRAIL annotate/INF1.merge-annotate.tsv | cut -f2-3,19,22 | tr ':' '\t' | sort -k2,2n -k3,3n
# grep TRAIL annotate/INF1.merge-annotate.tsv | cut -f2-3,19,22 | tr ':' '\t' | sort -k2,2n -k3,3n | \
# cut -f1 | grep -f - work/INF1.METAL | sort -k4,4n -k5,5n | cut -f1,2

Rscript -e '
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
  gz <- gzfile(file.path(INF,"METAL","TRAIL-1.tbl.gz"))
  TRAIL <- within(read.delim(gz,as.is=TRUE), {Z <- Effect/StdErr; P <- pvalue(Z); log10P <- -log10p(Z)}) %>%
           select(Chromosome,Position,MarkerName,Z,P,log10P)
  genes <- data.frame(chr=c("chr1", "chr3", "chr3", "chr11", "chr14", "chr17", "chr18", "chr19"),
                      snpid=c("chr1:196710916_C_T", "chr3:172274232_A_C", "chr3:186449122_A_G", "chr11:61549025_A_G",
                              "chr14:94844947_C_T", "chr17:64224775_C_T", "chr18:29804863_A_T", "chr19:44153100_A_G"),
                      snp=c("rs16840522", "rs574044675", "rs5030044", "rs174533", "rs28929474", "rs8178824", "rs654488", "rs4760"),
                      gene=c("CFHR1;CFHR3","TNFSF10","KNG1","MYRF","SERPINA1","APOH","MEP1B","PLAUR")

           )
  TRAIL <- left_join(TRAIL,genes,by=c("MarkerName"="snpid"),keep=TRUE) %>%
            mutate(MarkerName=ifelse(!is.na(snpid),gene,MarkerName)) %>%
            select(-c(chr,snpid,snp))
  save(TRAIL, genes, file=file.path(INF,"work","TRAIL.rda"))
  load(file.path(INF,"work","TRAIL.rda"))
  subset(TRAIL,!is.na(gene))
  log10p <- gap::log10p
  png("TRAIL-mhtplot.trunc.png", res=300, units="in", width=9, height=6)
  par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
  source(file.path(INF,"csd3","IL.12B-mhtplot.trunc.R"))
  mhtplot.trunc(TRAIL, chr="Chromosome", bp="Position", z="Z", snp="MarkerName",
                suggestiveline=FALSE, genomewideline=-log10(5e-10),
                cex.mtext=1.2, cex.text=1.2,
                annotatelog10P=-log10(5e-10), annotateTop = FALSE, highlight=with(genes,gene),
                mtext.line=3, y.brk1=30, y.brk2=70, delta=0.01, cex.axis=1.5, cex.x=2, cex=0.8, font=3, font.axis=1.5,
                y.ax.space=20,
                col = c("blue4", "skyblue")
  )
  dev.off()
'

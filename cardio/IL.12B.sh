#!/usr/bin/bash

export prot=IL.12B

function test()
{
  (
    gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
    head -1 | \
    cut -f1,2,10-12
    gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
    sed '1d' | \
    cut -f1,2,10-12 | \
    sort -k1,1n -k2,2n
  ) | \
  gzip -f > ${prot}.txt.gz

  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    library(gap)
    library(dplyr)
    chr_pos_b_se_p <- read.delim(paste0(prot,".txt.gz")) %>%
                      mutate(P=pvalue(Effect/StdErr)) %>%
                      select(Chromosome,Position,P)
    names(chr_pos_b_se_p) <- c("#CHROM","POS","P")
    head(chr_pos_b_se_p)
    f <- gzfile(file.path(INF,paste0(prot,".dat.gz")))
    write.table(chr_pos_b_se_p,file=f,quote=FALSE,row.names=FALSE)
  END

  (
    echo chromosome position
    grep ${prot} ${INF}/work/INF1.merge | \
    cut -f8,9 | \
    sort -k1,1n -k2,2n
  ) > ${prot}.annotate

  R --slave --vanilla --args \
    input_data_path=${INF}/${prot}.dat.gz \
    output_data_rootname=${prot} \
    custom_peak_annotation_file_path=${INF}/${prot}.annotate \
    reference_file_path=${INF}/cardio/turboman_hg19_reference_data.rda \
    pvalue_sign=5e-10 \
    plot_title="${protein} (${uniprot})" < ${INF}/cardio/turboman.r
}

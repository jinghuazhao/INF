#!/usr/bin/bash

export prot=IL.12B
export uniprot=P29460

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
  paste - <(echo SNP BHLHE40 LPP IL12B MHC SH2B3 FLT3 RAD51B TRAF3 | tr ' ' '\n') > ${INF}/work/${prot}.annotate

  (
    echo chromosome position SNP
    grep IL.12B ${INF}/annotate/INF1.merge-annotate.tsv | \
    cut -f3,22 | \
    sort -k1,1n -k2,2n | \
    sed 's/:/\t/'
  ) > ${INF}/work/${prot}.annotate

  export refgene_file_name=turboman_hg19_reference_data.rda
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    refgene_file_name <- Sys.getenv("refgene_file_name")
    load(file.path(INF,"cardio",refgene_file_name))
    library(dplyr)
    refgene_gene_coordinates_h19 <- refgene_gene_coordinates_h19 %>%
                                    mutate(gene_name=if_else(gene_name=="LOC285626","IL12B",gene_name)) %>%
                                    mutate(gene_name=if_else(gene_name=="PSORS1C3","MHC",gene_name)) %>%
                                    mutate(gene_name=if_else(gene_name=="SH2B3","ATXN2",gene_name)) %>%
                                    mutate(gene_name=if_else(gene_name=="FLT3","URAD",gene_name))
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

turboman

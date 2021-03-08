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
  ) > ${INF}/work/${prot}.annotate

  R --slave --vanilla --args \
    input_data_path=${INF}/work/${prot}.dat.gz \
    output_data_rootname=${prot} \
    custom_peak_annotation_file_path=${INF}/work/${prot}.annotate \
    reference_file_path=${INF}/cardio/turboman_hg19_reference_data.rda \
    pvalue_sign=5e-10 \
    plot_title="${prot} (${uniprot})" < ${INF}/cardio/turboman.r
}

turboman

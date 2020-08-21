#!/usr/bin/bash

export TMPDIR=/rds/user/jhz22/hpc-work/work
export dir=${INF}

R --no-save <<END
# https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf11
  libary(STRINGdb)
  STRINGdb$methods()
  STRINGdb$help("plot_network")
  INF1_merge_trans <- read.delim("work/INF1.merge.trans.vepoutput",skip=91,as.is=TRUE)
  string_db <- STRINGdb$new(version="11", species=9606)
  mapped <- string_db$map(INF1_merge_trans, "Gene", removeUnmappedRows=TRUE)
  hits <- with(mapped, STRING_id)
  pdf("string_db.pdf")
  string_db$plot_network(hits)
  enrichment <- string_db$get_enrichment(hits)
  annotations <- string_db$get_annotations(hits)
  clustersList <- string_db$get_clusters(hits)
  interactions <- string_db$get_interactions(hits)
  subnetwork <- string_db$get_subnetwork(hits)
  for(i in seq(1:3)){string_db$plot_network(clustersList[[i]])}
  dev.off()
  IL12B <- string_db$mp("IL12B")
  KITLG <- string_db$mp("KITLG")
  TNFSF10 <- string_db$mp("TNFSF10")
  string_db$get_neighbors(c(IL12B,KITLG,TNFSF10))
  string_db$get_interactions(c(IL12B,KITLG,TNFSF10))
END

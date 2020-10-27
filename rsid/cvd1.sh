#!/usr/bin/bash

export stables=https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-020-00287-2/MediaObjects/42255_2020_287_MOESM3_ESM.xlsx
R --no-save <<END
  stables <- Sys.getenv("stables")
  st1 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:23), rows=c(3:95))
  table(st1$INFI)
  write.table(names(table(st1$Short_annotation)),file="cvd1.txt",col.names=FALSE,quote=FALSE,row.names=FALSE)
  write.table(names(table(subset(st1,is.na(INFI))$Short_annotation)),file="cvd1-inf.txt",col.names=FALSE,quote=FALSE,row.names=FALSE)
  write.table(names(table(subset(st1,INFI=="Y")$Short_annotation)),file="inf.txt",col.names=FALSE,quote=FALSE,row.names=FALSE)
END

# https://zenodo.org/record/2615265#.X5f9zEdxeUk
export url=https://zenodo.org/record/2615265/files/
if [ ! -d ~/rds/results/public/proteomics/scallop-cvd1 ]; then mkdir ~/rds/results/public/proteomics/scallop-cvd1; fi
cat cvd1.txt | xargs -I {} bash -c "wget ${url}/{}.txt.gz -O ~/rds/results/public/proteomics/scallop-cvd1/{}.txt.gz"

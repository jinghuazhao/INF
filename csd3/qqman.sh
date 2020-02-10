# 10-2-2020 JHZ

R --no-save <<END
  require(gap.datasets)
  test <- mhtdata[c("chr","pos","p")]
  write.table(test,file="test.txt",row.names=FALSE,quote=FALSE)
  annotate <- unique(subset(mhtdata[c("chr","start","gene")],p<5e-8 & gene!=""))
  names(annotate) <- c("chromosome","position","nearest_gene_name")
  write.table(annotate,file="annotate.txt",row.names=FALSE,quote=FALSE)
END

R --slave --vanilla --args \
  input_data_path=test.txt \
  output_data_rootname=test_qq \
  plot_title="gap.datasets example" < turboqq.r

R --slave --vanilla --args \
  input_data_path=test.txt \
  output_data_rootname=test_man \
  custom_peak_annotation_file_path=annotate.txt \
  reference_file_path=turboman_hg19_reference_data.rda \
  pvalue_sign=5e-8 \
  plot_title="gap.datasets example" < turboman.r

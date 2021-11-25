# 27/3/2019 JHZ

# turbodata
R --no-save -q <<END
  library(gap.datasets)
  data <- with(mhtdata,cbind(chr,pos,p))
  glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3","PPP1R3B",
            "RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
  hdata <- subset(mhtdata,gene%in%glist)[c("chr","pos","p","gene")]
  names(hdata) <- c("chromosome","position","nearest_gene_name")
  gz <- gzfile("test.gz","w")
  write.table(data,gz,row.names=FALSE,quote=FALSE)
  close(gz)
  gz <- gzfile("test_annotation.gz","w")
  write.table(hdata,gz,row.names=FALSE,col.names=TRUE,quote=FALSE)
  close(gz)
END

# turboman
R --slave --vanilla --args \
  input_data_path=test.gz \
  output_data_rootname=test_man \
  custom_peak_annotation_file_path=test_annotation.gz \
  reference_file_path=turboman_hg19_reference_data.rda \
  pvalue_sign=5e-10 \
  plot_title="Manhattan plot" < turboman.r

# turboqq
R --slave --vanilla --args \
  input_data_path=test.gz \
  output_data_rootname=test_qq \
  plot_title="Q-Q plot" < turboqq.r

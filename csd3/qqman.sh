# 15-3-2020 JHZ

export turbo=/home/jhz22/cambridge-ceu/turbo

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
  plot_title="gap.datasets example" < ${turno}qq/turboqq.r

R --slave --vanilla --args \
  input_data_path=test.txt \
  output_data_rootname=test_man \
  custom_peak_annotation_file_path=annotate.txt \
  reference_file_path=${turbo}man/turboman_hg19_reference_data.rda \
  pvalue_sign=5e-8 \
  plot_title="gap.datasets example" < ${turbo}man/turboman.r

function CRP()
{
gunzip -c $SUMSTATS | cut -f1,3,8-11 | sed '1d' | grep -v NaN | awk '{split($1,a,":");print a[1],a[2],$6}' | grep -v X | sort -k1,1n -k2,2n > crp.dat

R --slave --vanilla --args \
  input_data_path=crp.dat \
  output_data_rootname=crp_qq \
  plot_title="CRP example" < ${turbo}qq/turboqq.r

R --slave --vanilla --args \
  input_data_path=crp.dat \
  output_data_rootname=crp_man \
  reference_file_path=${turbo}man/turboman_hg19_reference_data.rda \
  pvalue_sign=5e-8 \
  plot_title="CRP example" < ${turbo}man/turboman.r
}

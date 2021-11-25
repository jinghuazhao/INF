# 12-11-2019 JHZ

source tryggve/analysis.ini

echo "--> LD clumping"

export rt=$HOME/INF/sumstats/INTERVAL
ls $rt/*.gz | \
sed 's/INTERVAL.//g;s/.gz//g' | \
xargs -l basename | \
parallel -j4 --env rt -C' ' '
if [ -f work/{}.clumped ]; then rm work/{}.clumped; fi; \
gunzip -c $rt/INTERVAL.{}.gz | \
plink --bfile EUR \
      --clump $rt/INTERVAL.{}.gz \
      --clump-snp-field SNPID \
      --clump-field PVAL \
      --clump-kb 500 \
      --clump-p1 5e-10 \
      --clump-p2 0.0001 \
      --clump-r2 0 \
      --mac 50 \
      --out work/INTERVAL.{}; \
'
(
  grep CHR work/INTERVAL.*.clumped | \
  head -1
  grep -v CHR work/INTERVAL.*.clumped
) | \
sed 's|work/INTERVAL.||g;s/.clumped://g' | \
awk '(NF>1){$3="";print}' | \
awk '{$1=$1;if(NR==1)$1="prot";print}' > INTERVAL.clumped
R --no-save -q <<END
  require(gap)
  clumped <- read.table("INTERVAL.clumped",as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","BP","SNP","prot")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","uniprot")
  cistrans <- cis.vs.trans.classification(hits)
  sink("INTERVAL.clumped.out")
  with(cistrans,table)
  sink()
  sum(with(cistrans,table))
  pdf("INTERVAL.circlize.pdf")
  circos.cis.vs.trans.plot(hits="INTERVAL.clumped")
  dev.off()
END

echo "--> top signals"
# NOTE this was based on results from metal/20110325 without TRACKPOSITIONS

ls work/*clumped | \
sed 's|work/||g;s/.clumped//g' | \
xargs -l basename | \
parallel -j4 --env rt -C' ' '
( \
   grep -w {} st.bed | \
   awk -vOFS="\t" -vM=1000000 "{start=\$2-M;if(start<0) start=0;end=\$3+M;\$2=start;\$3=end};1" > st.tmp; \
   read chrom start end gene prot < st.tmp; \
   head -1 $rt/{}.clumped; \
   awk -vchr=$chrom "(NR > 1 && \$1==chr)" $rt/{}.clumped | \
   sort -k3,3 | \
   join -v1 -13 -21 - MHC.snpid | \
   sort -k2,2n -k3,3n \
) > $rt/{}.top
'

echo "--> GC lambda"

(
ls $rt/*.gz | \
sed 's/INTERVAL.//g;s/.gz//g' | \
xargs -l basename | \
parallel -j4 --env rt -C' ' '
  gunzip -c $rt/INTERVAL.{}.gz | \
  cut -f1,11 | \
  gzip -f > work/INTERVAL.{}.p.gz; \
  export protein={}; \
  R --no-save -q <<END
  library(gap); \
  rt <- Sys.getenv("rt"); \
  protein <- Sys.getenv("protein"); \
  gz <- gzfile(paste0("work/INTERVAL.",protein,".p.gz")); \
  p <- read.table(gz,as.is=TRUE,header=TRUE); \
  cat(protein,"GC.lambda=",gc.lambda(with(p,PVAL)),"\n")
END'
) > work/INTERVAL.lambda.log
grep GC.lambda work/INTERVAL.lambda.log | \
grep -v gc.lambda | \
sed 's/GC.lambda=//g' > work/INTERVAL.lambda.dat

echo "--> conditional analysis"

ls $rt/*.gz | \
sed 's/INTERVAL.//g;s/.gz//g' | \
xargs -l basename | \
parallel -j3 --env rt -C' ' '
( \
  echo SNP A1 A2 freq b se p N; \
  gunzip -c $rt/INTERVAL.{}.gz | \
  awk "(NR>1) { \
     snpid=\$1; \
     gsub(/chr/,\"\",snpid); \
     split(snpid,chrpos_a1_a2,\":\"); \
     chr=chrpos_a1_a2[1]; \
     split(chrpos_a1_a2[2],a,\"_\"); \
     pos=a[1]; \
     a1=toupper(\$6); \
     a2=toupper(\$7); \
     print \$1, a1, a2, \$8, \$9, \$10, \$11, \$5 \
  }" \
) > work/INTERVAL.{}.ma'

ls work/*.ma | \
sed 's/.ma//g' | \
xargs -l basename | \
sort | \
parallel -j3 --env rt -C' ' '
  if [ -f work/{}.jma.cojo ]; then rm work/{}.jma.cojo work/{}.ldr.cojo; fi; \
  gcta64 --bfile EUR --cojo-file work/{}.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
         --thread-num 3 --out work/{}
'
(
  grep SNP work/INTERVAL.*.jma.cojo | \
  head -1
  grep -v SNP work/INTERVAL.*.jma.cojo
) > work/INTERVAL.jma
(
  grep SNP work/INTERVAL.*.ldr.cojo | \
  head -1
  grep -v SNP work/INTERVAL.*.ldr.cojo
) > work/INTERVAL.ldr

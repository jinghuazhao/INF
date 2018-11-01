# 1-11-2018 JHZ

source analysis.ini

export rt=$HOME/INF/sumstats/INTERVAL

echo "--> top signals"

ls $rt/*.gz | \
xargs -l basename -s .gz | \
sed 's/INTERVAL.//g' | \
parallel -j4 --env rt -C' ' '
( \
  gunzip -c $rt/INTERVAL.{}.gz | \
  head -1; \
  gunzip -c $rt/INTERVAL.{}.gz | \
  awk "(NR > 1 && \$11 <= 5e-10 && \$8 > 0.0001)" | \
  sort -k1,1 | \
  join -v1 - MHC.snpid | \
  sort -k2,2n -k3,3n \
) > work/INTERVAL.{}.top
'
echo "--> LD clumping"

ls $rt/*.gz | \
xargs -l basename -s .gz | \
sed 's/INTERVAL.//g' | \
parallel -j4 --env rt -C' ' '
gunzip -c $rt/INTERVAL.{}.gz | \
awk "\$8 > 0.0001 && \$13" | \
cut -f1,11,13 | \
gzip -f > work/INTERVAL.{}.tmp; \
plink --bfile EUR1KG \
      --exclude MHC.snpid \
      --clump work/INTERVAL.{}.tmp \
      --clump-snp-field SNPID \
      --clump-field PVAL \
      --clump-kb 500 \
      --clump-p1 5e-10 \
      --clump-p2 0.0001 \
      --clump-r2 0.1 \
      --out work/INTERVAL.{}; \
rm work/INTERVAL.{}.tmp
'

rm -f work/INTERVAL.clumped
(
  grep CHR work/INTERVAL.*.clumped | \
  head -1
  grep -v CHR work/INTERVAL.*.clumped
) > work/INTERVAL.clumped

echo "--> GC lambda"

(
ls $rt/*.gz | \
xargs -l basename -s .gz | \
sed 's/INTERVAL.//g' | \
parallel -j4 --env rt -C' ' '
  gunzip -c $rt/INTERVAL.{}.gz | \
  cut -f1,11 | \
  gzip -f > work/INTERVAL.{}.p.gz; \
  export protein={}; \
  R --no-save -q <<END
  source("files/lambda.R"); \
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
xargs -l basename -s .gz | \
sed 's/INTERVAL.//g' | \
parallel -j3 --env rt -C' ' '
( \
echo SNP A1 A2 freq b se p N; \
gunzip -c $rt/INTERVAL.{}.gz | \
awk -vOFS="\t" "(NR>1) { \
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

ls work/*jma*cojo | \
xargs -l basename -s .jma.cojo | \
sort > INTERVAL.cojo.done

ls work/*.ma | \
xargs -l basename -s .ma | \
sort | \
join -v1 - INTERVAL.cojo.done | \
parallel -j3 --env rt -C' ' '
gcta64 --bfile EUR1KG --cojo-file work/{}.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
       --exclude-region-bp 6 30000000 5000 --thread-num 3 --out work/{}
'

rm -f work/INTERVAL.jma
(
  grep SNP work/INTERVAL.*.jma.cojo | \
  head -1
  grep -v SNP work/INTERVAL.*.jma.cojo
) > work/INTERVAL.jma

rm -f work/INTERVAL.ldr
(
  grep SNP work/INTERVAL.*.ldr.cojo | \
  head -1
  grep -v SNP work/INTERVAL.*.ldr.cojo
) > work/INTERVAL.ldr


echo "--> Q-Q, Manhattan, LocusZoom plots"

export p=IFN.gamma
ls $rt/INTERVAL.${p}.gz | \
xargs -l basename -s .gz | \
sed 's/INTERVAL.//g' | \
parallel -j1 --env rt -C' ' '
export protein={}; \
R --no-save <<END
protein <- Sys.getenv("protein");\
rt <- Sys.getenv("rt");\
gz <- gzfile(paste0(rt,"/INTERVAL.",protein,".gz"));\
qqman <- paste0("work/",protein,"-qqman.png");\
MarkerName <- "SNPID";\
PVAL <- "PVAL";\
source("files/qqman.R")
END'
(echo Chr Start End; echo 11 60739337 60786787) > st.bed
grep ${p} $HOME/INF/sumstats/INTERVAL.list | \
sed "s/INTERVAL_inf1_//g;s/_chr_merged.gz\*//g;s/___/ /g" | \
parallel -j${threads} -C' ' '
   gunzip -c /data/jampet/upload-20170920/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
   awk -f files/INTERVAL-lz.awk | \
   awk -f files/order.awk | \
   gzip -f > work/INTERVAL.{1}.gz'
'
awk 'NR>1' st.bed | \
parallel -j${threads} --env p -C' ' '
   gunzip -c work/INTERVAL.${p}.gz | \
   awk -vOFS="\t" -vchrom={1} -vStart={2} -vEnd={3} "(NR>1){ \
     snpid=\$1; \
     gsub(/chr/,\"\",snpid); \
     split(snpid,chrpos_a1_a2,\":\"); \
     chr=chrpos_a1_a2[1]; \
     split(chrpos_a1_a2[2],a,\"_\"); \
     pos=a[1]; \
     if (chr==chrom && pos >= Start && pos <= End) print \$1}" | \
   sort > st.tmp;
   gunzip -c work/INTERVAL.${p}.gz | \
   awk -vOFS="\t" "(NR>1 && \$11 != \"NA\" && \$15!=\".\") {print \$1,\$11,\$5,\$15}" | \
   sort -k1,1 | \
   join st.tmp - | \
   awk -vOFS="\t" "{if(NR==1) print \"MarkerName\", \"P-value\", \"Weight\";print \$4,\$2,\$3}"> work/INTERVAL.${p}.lz
'
awk 'NR>1' st.bed | \
parallel -j1 --env p -C' ' '
  cd work; \
  rm -f ld_cache.db; \
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal INTERVAL.${p}.lz \
            --plotonly --chr {1} --start {2} --end {3} --no-date --rundir .; \
  pdftopng chr{1}_{2}-{3}.pdf -r 300 INTERVAL.${p}; \
  xdg-open INTERVAL.${p}-000001.png; \
  cd -
'

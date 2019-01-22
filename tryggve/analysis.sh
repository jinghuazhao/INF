# 22-1-2019 JHZ

source tryggve/analysis.ini

echo "--> Q-Q/Manhattan/LocusZoom plots"

export rt=$HOME/INF
(
  echo -e "chrom\tstart\tend\tgene\tprot"
  sort -k2,2 $rt/inf1.list > inf1.tmp
  cut -f2,3,7-10 $rt/doc/olink.inf.panel.annot.tsv  | \
  awk -vFS="\t" -vOFS="\t" '(NR>1){
      gsub(/\"/,"",$0)
      if($2=="Q8NF90") $3="FGF5"
      if($2=="Q8WWJ7") $3="CD6"
      print
  }' | \
  sort -t$'\t' -k2,2 | \
  join -t$'\t' -j2 inf1.tmp - | \
  awk -vFS="\t" -vOFS="\t" '{print $5,$6,$7,$4,$2}' | \
  sort -k1,1n -k2,2n
) > st.bed
ls METAL/*-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j4 --env rt -C' ' 'export protein={}; R --no-save -q < $rt/tryggve/qqman.R'
ls METAL/*-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j3 -C' ' '
(
   echo -e "MarkerName\tP-value\tWeight"
   grep -w {} st.bed | \
   awk -vOFS="\t" -vM=1000000 "{start=\$2-M;if(start<0) start=0;end=\$3+M;\$2=start;\$3=end};1" > st.tmp; \
   read chrom start end gene prot < st.tmp; \
   gunzip -c METAL/{}-1.tbl.gz | \
   awk -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end "(\$1 == chr && \$2 >= start && \$2 <= end){split(\$3,a,\"_\");print a[1],\$12,\$14}"
)  > METAL/{}.lz'
ls METAL/*-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j1 -C' ' '
   grep -w {} st.bed | \
   awk -vOFS="\t" -vM=1000000 "{start=\$2-M;if(start<0) start=0;end=\$3+M;\$2=start;\$3=end};1" > st.tmp; \
   read chrom start end gene prot < st.tmp; \
   cd METAL; \
   rm -f ld_cache.db; \
   locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
             --plotonly --chr $chrom --start $start --end $end --no-date --rundir .; \
   mv chr${chrom}_${start}-${end}.pdf {}.lz.pdf; \
   pdftopng -r 300 {}.lz.pdf {}; \
   mv {}-000001.png {}.lz-1.png; \
   mv {}-000002.png {}.lz-2.png; \
   cd -
'
export IM_dir=/data/jinhua/ImageMagick-7.0.8-22/bin
$IM_dir/convert OPG.lz-1.png -resize 130% OPG.lz-3.png
$IM_dir/convert \( OPG.qq.png -append OPG.manhattan.png -append OPG.lz-3.png -append \) +append OPG-qml.png

echo "--> clumping"

export rt=$HOME/INF/METAL
ls METAL/*tbl.gz | \
sed 's/-1.tbl.gz//g' | \
xargs -l basename | \
parallel -j4 --env rt -C' ' '
if [ -f $rt/{}.clumped ]; then rm $rt/{}.clumped; fi; \
plink --bfile EUR \
      --clump $rt/{}-1.tbl.gz \
      --clump-snp-field MarkerName \
      --clump-field P-value \
      --clump-kb 500 \
      --clump-p1 5e-10 \
      --clump-p2 0.01 \
      --clump-r2 0 \
      --mac 50 \
      --out $rt/{}
'
(
  grep CHR $rt/*.clumped | \
  head -1
  grep -v CHR $rt/*.clumped
) | \
sed 's|'"$rt"'/||g;s/.clumped://g' | \
awk '(NF>1){$3="";print}' | \
awk '{$1=$1;if(NR==1)$1="prot";print}' > INF1.clumped
R --no-save -q <<END
  require(gap)
  clumped <- read.table("INF1.clumped",as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","BP","SNP","prot")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","uniprot")
  cistrans <- cis.vs.trans.classification(hits)
  sink("INF1.clumped.out")
  with(cistrans,table)
  sink()
  sum(with(cistrans,table))
  pdf("INF1.circlize.pdf")
  circos.cis.vs.trans.plot(hits="INF1.clumped")
  dev.off()
END

(
  gunzip -c METAL/4E.BP1-1.tbl.gz | \
  head -1
  awk 'NR>1' INF1.clumped | \
  cut -d' ' -f1,3 | \
  parallel -j1 -C' ' 'zgrep -H -w {2} METAL/{1}-1.tbl.gz'
) | \
sed 's|METAL/||g;s/-1.tbl.gz//g' > INF1.clumped.tbl

(
awk 'NR>1' INF1.clumped | \
cut -d' ' -f1,3 | \
parallel -j1 -C' ' '
  export direction=$(zgrep -w {2} METAL/{1}-1.tbl.gz | cut -f13)
  let j=1
  for i in $(grep "Input File" METAL/{1}-1.tbl.info | cut -d" " -f7)
  do
     export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
     if [ "$n" != "?" ]; then
        zgrep -H -w {2} $i
     fi
     let j=$j+1
  done
'
) | \
sed 's|/data/jinhua/INF/sumstats||g;s/.gz//g' > INF1.clumped.all

R -q --no-save <<END
  tbl <- read.delim("INF1.clumped.tbl",as.is=TRUE)
  tbl <- within(tbl, {
    prot <- sapply(strsplit(Chromosome,":"),"[",1)
    Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
  })
  all <- read.table("INF1.clumped.all",as.is=TRUE,col.names=c("SNPID", "CHR", "POS", "STRAND", "N",
         "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"))
  all <- within(all, {
    dir.study.prot <- sapply(strsplit(SNPID,":"),"[",1)
    p1 <- sapply(strsplit(SNPID,":"),"[",2)
    p2 <- sapply(strsplit(SNPID,":"),"[",3)
    MarkerName <- paste(p1,p2,sep=":")
    study <- sapply(strsplit(dir.study.prot,"/"),"[",2)
    study.prot <- sapply(strsplit(dir.study.prot,"/"),"[",3)
    substudy <- sapply(strsplit(study.prot,"[.]"),"[",1)
    p1 <- sapply(strsplit(study.prot,"[.]"),"[",2)
    p2 <- sapply(strsplit(study.prot,"[.]"),"[",3)
    prot <- ifelse(is.na(p2),p1,paste(p1,p2,sep="."))
  })
  require(rmeta)
  pdf("INF1.fp.pdf")
  for(i in 1:nrow(tbl)) if (i!=89 & i!=101)
  {
     p <- tbl[i,"prot"]
     m <- tbl[i,"MarkerName"]
     d <- gsub("[?]","",tbl[i,"Direction"])
     s <- unlist(strsplit(d,""))
     f <- as.numeric(paste0(s,1))
     A1 <- toupper(tbl[i,"Allele1"])
     A2 <- toupper(tbl[i,"Allele2"])
     print(paste0(i,"-",p,":",m))
     with(subset(all,prot==p & MarkerName==m), {
       xlim <- c(-1.5,1.5)
       e <- toupper(EFFECT_ALLELE)[1]
       r <- toupper(REFERENCE_ALLELE)[1]
       if (A1==e) {a1=A1;a2=A2;c=1;}
       else if(A1==r) {a1=A2;a2=A1;c=-1;}
       metaplot(BETA,SE,N,paste0(study," (",e,"/",r,")"),
                xlab="Effect distribution",ylab="",xlim=xlim,
                summn=tbl[i,"Effect"]*c,sumse=tbl[i,"StdErr"],sumnn=tbl[i,"N"])
       title(paste0(p," [",m," (",a1,"/",a2,")]"))
     })
  }
  dev.off()
END

echo "--> METAL results containing P-value=0"

awk '($5==0)' INF1.clumped | \
cut -d' ' -f1 | 
uniq > INF1.z
cat INF1.z | \
parall -j2 --env rt '
(
  export port={}
  gunzip -c $rt/{}-1.tbl.gz | \
  awk -vOFS="\t" "NR==1||\$12!=0";
  gunzip -c $rt/{}-1.tbl.gz | \
  awk -vOFS="\t" '(NR==1||\$12==0)' > {}.z; \
  R --no-save -q <<END
  prot <- Sys.getenv("prot");\
  metal <- read.delim(paste0(prot,".z"),as.is=TRUE);\
  library(Rmpfr) \
  metal <- within(metal,{P.value=format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE))});\
  write.table(metal,file=paste0(prot,".p"),sep=\"\\t\",row.names=FALSE,quote=FALSE);\
END
awk 'NR>1' {}.p;\
)'

echo "--> top signals"

export rt=$HOME/INF/METAL
ls METAL/*clumped | \
sed 's|METAL/||g;s/.clumped//g' | \
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

echo "--> COJO analysis"

export rt=$HOME/INF/METAL
ls METAL/*.tbl.gz | \
sed 's/-1.tbl.gz//g' | \
xargs -l basename | \
parallel -j3 --env rt -C' ' '
( \
  echo SNP A1 A2 freq b se p N; \
  gunzip -c $rt/{}-1.tbl.gz | \
  awk -vOFS="\t" "(NR>1 && \$14>50) {print \$3, \$4, \$5, \$6, \$10, \$11, \$12, \$14}" \
) > $rt/{}.ma
'
#1 Chromosome
#2 Position
#3 MarkerName
#4 Allele1
#5 Allele2
#6 Freq1
#7 FreqSE
#8 MinFreq
#9 MaxFreq
#10 Effect
#11 StdErr
#12 P-value
#13 Direction
#14 N

ls METAL/*.tbl.gz | \
sed 's/-1.tbl.gz//g' | \
xargs -l basename | \
parallel -j3 --env rt -C' ' '
if [ -f $rt/{}.jma.cojo ]; then rm $rt/{}.jma.cojo $rt/{}.ldr.cojo; fi; \
gcta64 --bfile EUR --cojo-file $rt/{}.ma --cojo-slct --cojo-p 5e-10 --cojo-collinear 0.01 --cojo-wind 500 \
       --maf 0.0001 --thread-num 3 --out $rt/{}
'
#1 MarkerName
#2 Allele1
#3 Allele2
#4 Freq1
#5 FreqSE
#6 MinFreq
#7 MaxFreq
#8 Effect
#9 StdErr
#10 P-value
#11 Direction
#12 CHR
#13 POS
#14 WEIGHT
(
  grep SNP $rt/*.jma.cojo | \
  head -1
  grep -v -w SNP $rt/*.jma.cojo
) | \
sed 's|'"$rt"'/||g;s/.jma.cojo:/\t/g' | \
awk -vOFS="\t" '{if(NR==1) $1="prot";print}' > INF1.jma

echo "--> clumping and cojo with LDetect approximately independent LD blocks"

awk 'NR>1{gsub(/chr/,"",$1);print}' tryggve/EURLD.bed > rlist-hg19
export rt=$HOME/INF
ls $rt/METAL/*tbl.gz | \
sed 's/-1.tbl.gz//g' | \
xargs -l basename | \
parallel -j6 --env rt -C' ' '
plink --bfile EUR \
      --clump $rt/METAL/{}-1.tbl.gz --clump-range rlist-hg19 \
      --clump-snp-field MarkerName \
      --clump-field P-value \
      --clump-p1 5e-10 --clump-p2 0.01 --clump-r2 0 \
      --mac 50 \
      --out $rt/LDBLOCK/{}'
awk '(NR>1){
  chr=$1;
  gsub(/chr/,"",chr);
  flanking=($3-$2)/2/1000
  centre=$2+flanking
  print sprintf("%d %d %d %s", chr, centre, flanking,$4);
}' tryggve/EURLD.bed > EURLD.region
for prot in $(ls $rt/METAL/*tbl.gz | sed 's/-1.tbl.gz//g' | xargs -l basename)
do
  export p=$prot
  # bruteforceclumpingbyregion
  cat EURLD.region | \
  parallel -j8 --env p --env rt -C' ' '
   gcta64 --bfile EUR --cojo-file $rt/METAL/$p.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
          --extract-region-bp {1} {2} {3} --thread-num 3 --out $rt/LDBLOCK/$p-{4}'
   (
     cat $rt/LDBLOCK/${p}*.jma.cojo | head -1
     awk "NR>1" $rt/LDBLOCK/${p}*.jma.cojo
   ) > $rt/${p}.jma
done

echo "--> contrast studies"

R --no-save -q <<END
z1z2 <- function(study1,study2,protein)
{
  file1 <- paste0("sumstats/",study1,"/",study1,".",protein,".gz")
  z1 <- gzfile(file1)
  interval <- read.table(z1,as.is=TRUE,header=TRUE)
  file2 <- paste0("sumstats/",study2,"/",study2,".",protein,".gz")
  z2 <- gzfile(file2)
  orcades <- read.table(z2,as.is=TRUE,header=TRUE)
  vars <- c("SNPID","BETA","SE")
  z <- merge(interval[vars],orcades[vars],by="SNPID")
  z <- within(z,{
    z1 <- BETA.x/SE.x
    z2 <- BETA.y/SE.y
  })
  with(z,cor(z1,z2))
  with(z,plot(z1,z2))
}
pdf("CXCL1.pdf")
z1z2("INTERVAL","ORCADES","IFN.gamma")
pdf()
END

# MAF x z plots
# BioConductor and CRAN packages

echo "--> Variant annotation"

export annovar_home=/services/tools/annovar/2018apr16
export humandb=$annovar_home/humandb
export example=$annovar_home/example
$annovar_home/annotate_variation.pl --geneanno -otherinfo -buildver hg19 $example/ex1.avinput $humandb/ --outfile ex1

# --- indirect features --- #

echo "--> remove duplicates"

# notes from https://www.biostars.org/p/264584/
# export LC_ALL=C 
#( \
#  grep '^#' input.vcf ; \
#  grep -v "^#" input.vcf | 
#  sort -t $'\t' -k1,1 -k2,2n -k4,4 | \
#  awk -F '\t' 'BEGIN {prev="";} {key=sprintf("%s\t%s\t%s",$1,$2,$4);if(key==prev) next;print;prev=key;}' \
#) 

function bruteforceclumpingbyregion()
{
  awk 'NR>1{gsub(/chr/,"",$1);print}' $rt/tryggve/EURLD.bed | \
  parallel -j8 --env p --env rt -C' ' '
  plink --bfile EUR \
      --chr {1} --from-bp {2} --to-bp {3} \
      --clump $rt/METAL/${p}-1.tbl.gz \
      --clump-snp-field MarkerName \
      --clump-field P-value \
      --clump-kb 500 \
      --clump-p1 5e-10 --clump-p2 0.01 --clump-r2 0 \
      --mac 50 \
      --out $rt/LDBLOCK/${p}-{4}'
  (
    cat $rt/LDBLOCK/${p}*.clumped | head -1
    awk "NR>1" $rt/LDBLOCK/${p}*.clumped
  ) > $rt/LDBLOCK/${p}.clumped
}

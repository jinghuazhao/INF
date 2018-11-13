# 12-11-2018 JHZ

source tryggve/analysis.ini

echo "--> Q-Q/Manhattan/LocusZoom plots"

ls METAL/*-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j2 --env rt -C' ' 'export protein={}; R --no-save -q < $rt/files/qqman.R'
# (echo Chr Start End; echo 4 73649784 76033785) > st.bed
(
  echo -e "chrom\tstart\tend\tgene\tprot"
  sort -k2,2 inf1.list > inf1.tmp
  cut -f2,3,7-10 doc/olink.inf.panel.annot.tsv  | \
  awk -vOFS="\t" '(NR>1){
      gsub(/\"/,"",$0)
      if($2=="Q8NF90") $3="FGF5"
      if($2=="Q8WWJ7") $3="CD6"
      print
  }' | \
  sort -k2,2 | \
  join -j2 inf1.tmp - | \
  awk -vOFS="\t" '{print $5,$6,$7,$4,$2}'
) > st.bed
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
   mv chr${chrom}_${start}-${end}.pdf {}.pdf; \
   pdftopng -r 300 {}.pdf {}; \
   cd -
'
# convert \( OPG-qqman-000001.png -append OPG-qqman-000002.png -append OPG-000001.png -append \) +append OPG.png

echo "--> 1000Genomes reference data"

# individual genotypes
export EUR=/services/tools/locuszoom/1.4/data/1000G/genotypes/2014-10-14/EUR

seq 22 | \
parallel 'plink --bfile $EUR/chr{} --recode vcf bgz --out chr{}'
seq 22 | \
parallel 'tabix -f -p vcf chr{}.vcf.gz'
seq 22 | \
awk -vp=$PWD '{print p "/chr" $1 ".vcf.gz"}' > EUR.list
vcf-concat --files EUR.list | \
bcftools norm -d both - | \
bcftools convert - -O z -o EUR.vcf.gz
plink --vcf EUR.vcf.gz --make-bed --out EUR
awk -vFS="\t" '
{
  CHR=$1
  POS=$4
  a1=$5
  a2=$6
  if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
  else snpid="chr" CHR ":" POS "_" a1 "_" a2
  print snpid, $2

# 1 CHR
# 2 SNPID
# 3 0
# 4 POS
# 5 A1
# 6 A2

}' EUR.bim > EUR.snpid

plink --bfile EUR --update-name EUR.snpid 1 2 --make-bed --out EUR1KG
plink --bfile EUR1KG --chr 6 --from-mb 25 --to-mb 35 --make-bed --out MHC
cut -f2 MHC.bim > MHC.snpid

echo "--> clumping"

export rt=$HOME/INF/METAL

ls METAL/*tbl.gz | \
sed 's/-1.tbl.gz//g' | \
xargs -l basename | \
parallel -j4 --env rt -C' ' '
plink --bfile EUR1KG \
      --exclude MHC.snpid \
      --clump $rt/{}-1.tbl.gz \
      --clump-snp-field MarkerName \
      --clump-field P-value \
      --clump-kb 500 \
      --clump-p1 5e-10 \
      --clump-p2 0.0001 \
      --clump-r2 0.1 \
      --out $rt/{}
'

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

echo "--> COJO analysis, --cojo-wind 10000"

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
gcta64 --bfile EUR1KG --cojo-file $rt/{}.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
       --exclude-region-bp 6 30000000 5000 --thread-num 3 --out $rt/{}
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

echo "--> clumping and cojo with LDetect approximately independent LD blocks"

awk '(NR>1){
  chr=$1;
  gsub(/chr/,"",chr);
  flanking=($3-$2)/2/1000
  centre=$2+flanking
  print sprintf("%d %d %d %s", chr, centre, flanking,$4);
}' tryggve/EURLD.bed > EURLD.region
export rt=$HOME/INF/METAL
for p in $(ls METAL/*tbl.gz | sed 's/-1.tbl.gz//g' | xargs -l basename)
do
  awk 'NR>1' tryggve/EURLD.bed | \
  parallel --env p --env rt -C' ' '
  ( \
    plink --bfile EUR1KG \
      --chr {1} \
      --from-bp {2} \
      --to-bp {3} \
      --exclude MHC.snpid \
      --clump $rt/{}-1.tbl.gz \
      --clump-snp-field MarkerName \
      --clump-field P-value \
      --clump-kb 500 \
      --clump-p1 5e-10 \
      --clump-p2 0.0001 \
      --clump-r2 0.1 \
      --out $rt/${p}-{4}; \
    if [ -f ${p}-{4}.clumped ]; then awk "NR>1" $rt/$p-{4}.clumped; fi; \
   ) > $rt/${p}.clumped; \
  cat EURLD.region | \
  parallel --env p --env rt -C' ' '
   ( \
     gcta64 --bfile KORA2 --cojo-file $rt/$p.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
            --extract-region-bp {1} {2} {3} --exclude-region-bp 6 30000000 5000 --thread-num 3 --out METAL/$p-{4}; \
     if [ -f ${p}-{4}.jma.cojo ]; then awk "NR>1" $rt/$p-{4}.jma.cojo; fi; \
   ) > $rt/${p}.jma
  '
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

echo "--> SNP/gene databases and indices"

export hg19=/services/tools/locuszoom/1.4/data/database/locuszoom_hg19.db
sqlite3 $hg19 <<END
-- locuszoom 1.4 tables and indices
.tables
.indices
.separator "\t"
.header on
.output gencode.tsv
select * from gencode;
.output recomb_rate.tsv
select * from recomb_rate;
.output refFlat.tsv
select * from refFlat;
.output refsnp_trans.tsv
select * from refsnp_trans;
.output snp_pos.tsv
select * from snp_pos;
.output snp_set.tsv
select * from snp_set;
END

gzip -f snp_pos.tsv

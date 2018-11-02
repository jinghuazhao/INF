# 31-10-2018 JHZ

source analysis.ini

echo "--> Q-Q/Manhattan/LocusZoom plots"

export p=IFN.gamma
ls METAL/${p}-1.tbl.gz | \
sed 's|METAL/||g;s/-1.tbl.gz//g' | \
parallel -j$threads -C' ' '
export protein={}; \
R --no-save <<END
protein <- Sys.getenv("protein");\
gz <- gzfile(paste0("METAL/",protein,"-1.tbl.gz"));\
qqman <- paste0("METAL/",protein,"-qqman.png");\
MarkerName <- "MarkerName";\
PVAL <- "P.Value";\
source("files/qqman.R");\
END'
(echo Chr Start End; echo 4 73649784 76033785) > st.bed
awk 'NR>1' st.bed | \
parallel -j${threads} --env p -C' ' '
   gunzip -c METAL/${p}-1.tbl.gz | \
   awk -vOFS="\t" -vchrom={1} -vStart={2} -vEnd={3} "(NR>1){
     snpid=\$1; \
     gsub(/chr/,\"\",snpid); \
     split(snpid,chrpos_a1_a2,\":\"); \
     chr=chrpos_a1_a2[1]; \
     split(chrpos_a1_a2[2],a,\"_\"); \
     pos=a[1]; \
     if (chr == chrom && pos >= Start && pos <= End) print \$1}" | \
   sort > st.tmp;
   gunzip -c METAL/${p}-1.tbl.gz | \
   awk -vOFS="\t" "(NR>1) {print \$1,\$10,\$14}" | \
   sort -k1,1 | \
   join st.tmp - | \
   awk -vOFS="\t" "{if(NR==1) print \"MarkerName\", \"P-value\", \"Weight\";print \$1,\$2,\$3}"> METAL/${f}.lz
'
awk 'NR>1' st.bed | parallel -j1 --env p -C' ' '
  rm -f ld_cache.db; \
  export f=chr{1}_{2}_{3}; \
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal METAL/$f.lz \
            --plotonly --chr {1} --start {2} --end {3} --no-date --rundir .; \
  pdftopng chr{1}_{2}-{3}.pdf -r 300 $p; \
  xdg-open ${p}-000001.png
'

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
xargs -l basename -s -1.tbl.gz | \
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

echo "--> COJO analysis, --cojo-wind 10000"

export rt=$HOME/INF/METAL
ls METAL/*.tbl.gz | \
xargs -l basename -s -1.tbl.gz | \
parallel -j3 --env rt -C' ' '
( \
echo SNP A1 A2 freq b se p N; \
gunzip -c $rt/{}-1.tbl.gz | \
awk -vOFS="\t" "(NR>1 && \$14>0) { \
   snpid=\$1; \
   gsub(/chr/,\"\",snpid); \
   split(snpid,chrpos_a1_a2,\":\"); \
   chr=chrpos_a1_a2[1]; \
   split(chrpos_a1_a2[2],a,\"_\"); \
   pos=a[1]; \
   a1=toupper(\$2); \
   a2=toupper(\$3); \
   print \$1, a1, a2, \$4, \$8, \$9, \$10, \$14, chr, pos \
}" | \
sort -k9,9n -k10,10n | \
cut -f1-8 --output-delimiter=" " \
) > $rt/{}.ma; \

ls METAL/*.tbl.gz | \
xargs -l basename -s -1.tbl.gz | \
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

echo "--> LDetect, approximate LD blocks"

for p in $(ls METAL/*tbl.gz | xargs -l basename -s -1.tbl.gz)
do
awk '(NR>1){
  chr=$1;
  gsub(/chr/,"",chr);
  region=sprintf("%d %d %s", chr, $2+($3-$2)/2, $4);
  print region
}' EURLD.bed | \
parallel --dry-run --env rt --env p -C' ' '
  rm -rf work/$p.jma; \
  touch work/$p.jma; \
  gcta64 --bfile KORA2 --cojo-file $rt/$p.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
         --extract-region-bp {1} {2} 1 --exclude-region-bp 6 30000000 5000 --thread-num 3 --out work/$p-{3}; \
  awk "NR>1" work/$p-{3}.jma.cojo >> work/$p.jma'
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

echo "--> add positions to METAL .tbl files"

gunzip -c snp_pos.tsv.gz | \
awk 'NR>1' | \
sort -k1,1 > snp_pos.noh

cut -f1 inf1.list | \
parallel -j6 -C' ' 'awk -vOFS="\t" "BEGIN{print \"SNP\",\"CHR\",\"BP\",\"P\",\"zscore\"}" > METAL/{}.tsv'
ls METAL/*.tbl | \
xargs -l basename -s -1.tbl | \
parallel -j6 -C' ' '
awk "NR>1" METAL/{}-1.tbl | \
sort -k1,1 | \
join -j1 snp_pos.noh - | \
sort -k2,2n -k3,3n | \
awk -vOFS="\t" "{print \$1, \$2, \$3, \$8, \$7}" >> METAL/{}.tsv'

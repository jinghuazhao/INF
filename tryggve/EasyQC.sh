# 3-12-2018 JHZ

source tryggve/analysis.ini

# 1000Genomes phase 3 from EasyQC
cd /data/jinhua/data/EasyQC
gunzip -c 1000GP_p3v5_legends_rbind.noDup.noMono.noCnv.noCnAll.afref.EUR.txt.gz | \
awk -vOFS="\t" '
{
  if(NR==1) print "SNP","CHR","POS","MINOR","MAJOR","MAF";
  else
  {
     split($1,a,":");
     a1=$2
     a2=$3
     if (a1>a2) snpid="chr" $1 "_" a2 "_" a1;
     else snpid="chr" $1 "_" a1 "_" a2
     if ($4<1-$4) {minor=$2;major=$3;maf=$4}
     else {minor=$3;major=$2;maf=1-$4}
     print snpid,a[2],a[1],minor,major,maf
  }
}' | \
gzip -f > 1KGp3v5.tsv.gz

R --no-save -q <<END
   z <- gzfile("1KGp3v5.tsv.gz")
   allele_ref_std <- read.table(z,header=TRUE,as.is=TRUE)
   save(allele_ref_std,file="1KGp3v5.RData")
END
cd -

# cptid   ea      oa      eaf
# 10:100000122    T       A       0.997018

#> load("/data/jinhua/data/QCGWAS/hapmap.RData")
#> ls()
#[1] "allele_ref_std"
#> head(allele_ref_std)
#         SNP CHR    POS MINOR MAJOR MAF SOURCE               DATE_ADDED
#1 rs10399749   1  45162     T     C   0 hapmap Thu Nov 22 16:08:25 2018
#2  rs4030303   1  72434     A     G   0 hapmap Thu Nov 22 16:08:25 2018
#> dim(allele_ref_std)
#[1] 4026413       8

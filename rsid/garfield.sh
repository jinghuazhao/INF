#!/usr/bin/bash

wget -qO- https://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield-data.tar.gz | \
tar xfz -

# HLA-treated SNPs
ls sentinels/*.p | \
xargs -l sed '1d;s/chr//' | \
sort -k1,1n -k3,3n | \
cut -f2 --complement > work/garfield.dat

# all significant SNPs
(
  ls METAL/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz* | \
  parallel -j3 -C' ' 'zcat METAL/{}-1.tbl.gz | awk " NR>1 && \$12<-5"'
) | sort -k1,1n -k2,2n > work/garfield.dat

# garfield-create-input-gwas.sh
# the column in GWAS file containing chormosome information
# the column in GWAS file containing genomic position information
# the column in GWAS file containing GWAS p-value information
# name of directory for GWAS trait to be created 
# name of file containing GWAS summary statistics
chrcol=1
poscol=2
pvalcol=12
TRAITNAME=INF1
GWASFILENAME=$INF/work/garfield.dat

# output directory to be used as input for GARFIELD analysis
OUTDIR=garfield-data/pval/$TRAITNAME
mkdir -p $OUTDIR

for CHR in {1..22}
do
  echo $CHR
  awk -v chr=$CHR -v chrcol=$chrcol -v poscol=$poscol -v pvalcol=$pvalcol '
      $chrcol==chr {print $poscol,10^$pvalcol}' $GWASFILENAME | \
  sort -k1n > $OUTDIR/chr$CHR
done

module load gcc/6

R --no-save <<END
  library(garfield)
  garfield.run("INF1", data.dir="garfield-data", trait="INF1", run.option = "prep", chrs = 1:22)
  n.perm <- 100000
  e <- c(5:10,100)
  garfield.run("INF1", data.dir="garfield-data", run.option = "perm", nperm = n.perm,
               thresh = 10^-e, pt_thresh = 10^-e,  maf.bins = 5, tags.bins = 5, tss.bins = 5,
               prep.file = "INF1.prep", optim_mode = TRUE, minit = 100, thresh_perm = 0.0001)
  garfield.plot("INF1.perm", num_perm = n.perm,
                output_prefix = "INF1", plot_title = "SCALLOP/INF1", filter = 10, tr = Inf)
  p <- read.table("INF1.perm",header=TRUE,as.is=TRUE)
  dim(p)
  attach(p)
  table(Tissue)
  length(table(Tissue))
  table(Type)
  length(table(Type))
  table(Celltype)
  length(table(Celltype))
  table(Category)
  length(table(Category))
  detach(p)
END

cd /rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/R/garfield-v2
garfield-INF1

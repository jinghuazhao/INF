#!/usr/bin/bash

function download()
{
  wget -qO- https://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield-data.tar.gz | \
  tar xfz -
}

function HLA()
# HLA-treated SNPs
{
  ls sentinels/*.p | \
  xargs -l sed '1d;s/chr//' | \
  sort -k1,1n -k3,3n | \
  cut -f2 --complement > work/garfield.dat
}

function all_snps()
{
# all significant SNPs
(
  ls ${INF}/METAL/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz* | \
  parallel -j3 -C' ' 'zcat METAL/{}-1.tbl.gz | awk " NR>1 && \$12<-5"'
) | sort -k1,1n -k2,2n > ${INF}/work/garfield.dat
}

function garfield_input()
{
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
}

function v1()
# v1
{
module load gcc/6
R --no-save <<END
  library(garfield)
  garfield.run("INF1", data.dir="garfield-data", trait="INF1", run.option = "prep", chrs = 1:22)
  n.perm <- 1e50
  e <- 50
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
}

function v2()
# v2
{
  cd ${HPC_WORK}/garfield-v2
  garfield-INF1
}

# Celltypes
cut -d' ' -f15 garfield.test.INF1.out | sort | grep -v -w -e Celltype -e NA | uniq | wc -l

# --- by proteins

function protein_snps()
{
# all significant SNPs
  chmod +x ${INF}/METAL/*-1.tbl.gz
  ls ${INF}/METAL/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz* | \
  parallel -j10 --env INF -C' ' 'zcat ${INF}/METAL/{}-1.tbl.gz | awk "NR>1 && \$12<-5" | \
                                 sort -k1,1n -k2,2n > ${INF}/garfield/garfield-{}.dat'
}

function protein_input()
{
  ls ${INF}/METAL/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz* | parallel -j10 --env INF -C' ' '
  export chrcol=1
  export poscol=2
  export pvalcol=12
  export TRAITNAME={}
  export GWASFILENAME=${INF}/garfield/garfield-{}.dat

# output directory to be used as input for GARFIELD analysis
  export OUTDIR=${INF}/garfield-data/pval/$TRAITNAME
  mkdir -p ${OUTDIR}

  for CHR in {1..22}
  do
    echo {} - $CHR
    awk "\$chrcol==chr {print \$poscol,10^\$pvalcol}" chr=$CHR chrcol=$chrcol poscol=$poscol pvalcol=$pvalcol $GWASFILENAME | \
    sort -k1n > ${OUTDIR}/chr$CHR
  done
  '
}

function protein_v2()
{
  cd ${HPC_WORK}/garfield-v2
  garfield-prot
}

function extract_v2()
# ID PThresh OR Pvalue Beta SE CI95_lower CI95_upper NAnnotThesh NAnnot NThresh N linkID Annotation Celltype Tissue Type Category
{
(
  ls ${INF}/METAL/*-1.tbl.gz | \
  xargs -l basename -s -1.tbl.gz* | parallel -j10 --env INF -C' ' '
  export prot={}
  export SCALLOP=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/
  if [ -f $SCALLOP/R/garfield-data/output/{}/garfield.test.{}.out ]; then
R --no-save -q <<END
  options(echo=FALSE,width=350)
  SCALLOP <- Sys.getenv("SCALLOP")
  prot <- Sys.getenv("prot")
  f <- paste0("garfield.test.",prot,".out")
  test <- read.table(file.path(SCALLOP,"R","garfield-data","output",prot,f),header=TRUE)
  sset <- subset(test,!is.na(Pvalue)&Pvalue<1e-3)[,c(1:8,14:18)]
  if (nrow(sset)>0) print(data.frame(prot,sset),row.names=FALSE)
END
  fi
  '
) | grep -v options | awk 'NR==1 || (NR>1 && !/OR/)' > ${INF}/garfield/garfield-3.txt

R --no-save -q <<END
  testfun <- function()
  {
  # clear but less elegant.
    par(mfrow=c(1,2))
    with(INF1,gap::ESplot(INF1[c("ID","logOR","CI95_lower","CI95_upper")],SE=FALSE,logscale=FALSE,xlim=c(0.1,0.9),v=0))
    title("All protein effects")
    with(prot3,gap::ESplot(prot3[c("ID","logOR","CI95_lower","CI95_upper")],SE=FALSE,logscale=FALSE,xlim=c(-0.5,5.5),v=10))
    title("Protein specific effects")
  }
#
  library(ggplot2)
  esplot <- function(data,sep,xlim,breaks,title)
  {
  # Somehow there are extra empty lines.
    p <- ggplot(data=data, aes(y=index, x=Beta, xmin=CI95_lower, xmax=CI95_upper))+
    geom_point()+
    geom_errorbarh()+
    scale_x_continuous(limits = xlim, breaks = breaks, name=expression(""))+
    scale_y_continuous(name = "", breaks=with(data,index), labels = with(data,ID), trans="reverse")+
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
    facet_grid(Category~., scales= "free", space="free")+
    ggtitle(title)+
    theme_minimal()+
    theme(text=element_text(family="Times",size=18, color="black"))+
    theme(panel.spacing = unit(1, "lines"))
    p
  }
  INF <- Sys.getenv("INF")
  library(dplyr)
  library(ggforestplot)
  library(ggplot2)
  INF1 <- read.table(file.path(INF,"garfield-data","output","INF1","garfield.test.INF1.out"),header=TRUE) %>%
          filter(Pvalue<=1e-6 & !is.na(Tissue)) %>%
          mutate(ID=paste(Tissue,Celltype,sep=":",if_else(Category=="Hotspots","HS","HM")),logOR=log(OR),index=1:n(),prot="") %>%
          arrange(desc(Category))
  prot3 <- read.table(file.path(INF,"garfield","garfield-3.txt"),header=TRUE) %>%
           mutate(prot=if_else(prot=="FGF.5","FGF-5",prot),
                  ID=paste(prot,sep=":",Tissue,Celltype,if_else(Category=="Hotspots","HS","HM")),logOR=log(OR),index=1:n())
# p1 <- esplot(INF1,sep="",xlim=c(-0.1,1.5),breaks=seq(-0.1,0.85,0.2),title="All protein effects")
# p2 <- esplot(prot3,sep="-",xlim=c(-0.5,6),breaks=seq(-0.5,5.5,by=2),title="Protein-specific effects")
#
  p1 <- ggforestplot::forestplot(INF1, name = ID, estimate = Beta, se = SE)
  p2 <- ggforestplot::forestplot(prot3, name = ID, estimate = Beta, se = SE)
  ggsave(p1,filename=file.path(INF,"garfield","garfield-INF1.png"),device="png")
  ggsave(p2,filename=file.path(INF,"garfield","garfield-prot.png"),device="png")
END
}

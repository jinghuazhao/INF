# 21-9-2019 JHZ

export TMPDIR=/rds/user/jhz22/hpc-work/work
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export tag=_nold

for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do

echo $p
export p=${p}
R --no-save -q <<END
  require(Rmpfr)
  INF <- Sys.getenv("INF")
  p <- Sys.getenv("p")
  tag <- Sys.getenv("tag")
  d <- read.delim(paste0(INF,"/sentinels/",p,tag,".p"))
  d <- within(d, {
    CHR <- sub("chr","",Chrom)
    POS <- End
    EPACTS <- sprintf("%s:%d_%s/%s",CHR,POS,toupper(Allele1),toupper(Allele2))
    PVALUE <- as.numeric(2*pnorm(mpfr(-abs(Effect/StdErr),100),lower.tail=TRUE,log.p=FALSE))
    log10p <- log.P./log(10)
    prot <- p
  })[c("MarkerName","EPACTS","CHR","POS","PVALUE","Effect","StdErr","log10p","prot")]
  write.table(d,file=paste0(INF,"/work/",p,".p"),quote=FALSE,row.names=FALSE,sep="\t")
END

function dummy()
{
awk -v OFS="\t" -v prot=$p '
{
   if (NR==1) print "MarkerName","EPACTS","CHR","POS","PVALUE","Effect","StdErr","log10p","prot";
   else {
     gsub(/chr/,"",$1)
     EPACTS=$1 ":" $3 "_" toupper($5) "/" toupper($6)
     print $4,EPACTS,$1,$3,exp($13),$11,$12,$13/log(10),prot
   }
}' ${INF}/sentinels/${p}${tag}.p > work/${p}.p
}

swiss --assoc ${INF}/work/${p}.p \
      --variant-col EPACTS \
      --chrom-col CHR \
      --pos-col POS \
      --pval-col PVALUE \
      --include-cols MarkerName Effect StdErr log10P --trait ${p} --skip-gwas \
      --dist-clump --clump-dist 1000000 --clump-p 5e-10 --out ${INF}/work/${p}
done
(
  cat ${INF}/work/*clump | head -1
  for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do awk 'NR>1' ${INF}/work/${p}.clump; done
) > INF1.hits

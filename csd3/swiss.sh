# 22-9-2019 JHZ

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
    BP <- End
    SNP <- MarkerName
    EPACTS <- sprintf("%s:%d_%s/%s",CHR,POS,toupper(Allele1),toupper(Allele2))
    PVALUE <- as.numeric(2*pnorm(mpfr(-abs(Effect/StdErr),100),lower.tail=TRUE,log.p=FALSE))
    log10p <- log.P./log(10)
    prot <- p
  })[c("MarkerName","EPACTS","CHR","POS","PVALUE","Effect","StdErr","log10p","prot","BP","SNP")]
  write.table(d,file=paste0(INF,"/work/",p,".p"),quote=FALSE,row.names=FALSE,sep="\t")
END

swiss --assoc ${INF}/work/${p}.p \
      --variant-col EPACTS \
      --chrom-col CHR \
      --pos-col POS \
      --pval-col PVALUE \
      --include-cols MarkerName Effect StdErr log10P BP SNP --trait ${p} --skip-gwas \
      --dist-clump --clump-dist 1000000 --clump-p 5e-10 --out ${INF}/work/${p}
done
(
  cat ${INF}/work/*clump | head -1
  for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do awk 'NR>1' ${INF}/work/${p}.clump; done
) > INF1.hits

R --no-save -q <<END
  require(gap)
  tag <- Sys.getenv("tag")
  rt <- paste0("INF1")
  clumped <- read.table(paste0(rt,".hits"),as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","POS","MarkerName","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
  cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
  cis.vs.trans <- with(cistrans,data)
  write.table(cis.vs.trans,file=paste0(rt,".hits.cis.vs.trans"),row.names=FALSE,quote=TRUE)
  cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
  write.table(cis,file=paste0(rt,".hits.cis"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  sink(paste0(rt,".hits.out"))
  with(cistrans,table)
  sink()
  with(cistrans,total)
  pdf(paste0(rt,".hits.circlize.pdf"))
  circos.cis.vs.trans.plot(hits=paste0(rt,".hits"),inf1,"uniprot")
  dev.off()
END

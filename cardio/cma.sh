# 2-7-2019 JHZ

(
  cat dist/cojo/*cma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  cut -d' ' -f1,4 dist/INF1_nold.sentinels | \
  parallel -j1 -C' ' '
    if [ -f dist/cojo/{1}-{2}.cma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "\$NF <= 5e-10 {print prot, snpid, \$0}" dist/cojo/{1}-{2}.cma.cojo
    fi'
) > INF1.cma

(
  awk -vOFS="," 'BEGIN{print "prot","CHR","BP","SNP","l","u","d","log10p","Groupid", "Type"}'
  (
  R --no-save -q <<\ \ END
  options(echo=FALSE)
  require(gap)
  require(reshape)
  cma <- read.delim("INF1.cma", as.is=TRUE)
  cma <- rename(cma,c(Chr="Chrom",bp="End",bC="Effect",bC_se="StdErr",SNP="MarkerName"))
  tag <- "_cma"

  for (protein in with(cma, unique(prot)))
  {
    p <- subset(cma, prot==protein)
    chrs <- with(p,unique(Chrom))
    for(chr in chrs)
    {
      ps <- subset(p,Chrom==chr)
      row.names(ps) <- 1:nrow(ps)
      sentinels(ps, protein, 1)
    }
  }
  END
  ) | \
  awk -vFS="," -vOFS="," '!/option/{
       SNPID=$2
       split(SNPID,a,":")
       split(a[2],b,"_")
       gsub(/chr/,"",a[1])
       $1=$1 "," a[1] "," b[1]
       print
  }'
) | \
sed 's/,/\t/g' > INF1.cma.sentinels

R --no-save -q <<END
    require(gap)
    clumped <- read.delim("INF1.cma.sentinels",as.is=TRUE,header=TRUE)
    hits <- merge(clumped[c("CHR","BP","SNP","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
    cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
    cis.vs.trans <- with(cistrans,data)
    write.table(cis.vs.trans,file="INF1.cma.sentinels.cis.vs.trans",row.names=FALSE,quote=TRUE)
    cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
    write.table(cis,file="INF1.cma.sentinels.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
    sink("INF1.cma.sentinels.out")
    with(cistrans,table)
    sink()
    with(cistrans,total)
    pdf("INF1.cma.sentinels.circlize.pdf")
    circos.cis.vs.trans.plot(hits="INF1.cma.sentinels",inf1,"uniprot")
    dev.off()
END


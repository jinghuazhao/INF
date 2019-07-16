# 16-7-2019 JHZ

(
  cat INTERVAL/*cma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  cut -d' ' -f1,4 INTERVAL/INTERVAL_nold.sentinels | \
  parallel -j1 -C' ' '
    if [ -f INTERVAL/{1}-{2}.cma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "\$NF <= 5e-10 {print prot, snpid, \$0}" INTERVAL/{1}-{2}.cma.cojo
    fi'
) > INTERVAL.cma

(
  awk -vOFS="\t" 'BEGIN{print "prot","CHR","BP","SNP","l","u","d","log10p","Groupid", "Type"}'
  (
  R --no-save -q <<\ \ END
  options(echo=FALSE)
  require(gap)
  require(reshape)
  cma <- read.delim("INTERVAL.cma", as.is=TRUE)[c("prot","Chr","bp","bC","bC_se","SNP")]
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
      sentinels(ps, protein, 1, sep="\t")
    }
  }
  END
  ) | \
  awk -vFS="\t" -vOFS="\t" '!/option/{
       SNPID=$2
       split(SNPID,a,":")
       split(a[2],b,"_")
       gsub(/chr/,"",a[1])
       $1=$1 "\t" a[1] "\t" b[1]
       print
  }'
) > INTERVAL.cma.sentinels

R --no-save -q <<END
    require(gap)
    clumped <- read.delim("INTERVAL.cma.sentinels",as.is=TRUE,header=TRUE)
    hits <- merge(clumped[c("CHR","BP","SNP","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
    cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
    cis.vs.trans <- with(cistrans,data)
    write.table(cis.vs.trans,file="INTERVAL.cma.sentinels.cis.vs.trans",row.names=FALSE,quote=TRUE)
    cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
    write.table(cis,file="INTERVAL.cma.sentinels.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
    sink("INTERVAL.cma.sentinels.out")
    with(cistrans,table)
    sink()
    with(cistrans,total)
    pdf("INTERVAL.cma.sentinels.circlize.pdf")
    circos.cis.vs.trans.plot(hits="INTERVAL.cma.sentinels",inf1,"uniprot")
    dev.off()
END

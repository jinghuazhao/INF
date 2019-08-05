# 5-8-2019 JHZ
(
  cat work/*jma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  cut -d' ' -f1,4 work/INF1_nold.sentinels | \
  parallel -j1 -C' ' '
    if [ -f work/{1}-{2}.jma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "NR > 1 {print prot, snpid, \$0}" work/{1}-{2}.jma.cojo
    fi'
) > INF1.jma
sed 's/Chr/CHR/g;s/bp/BP/g' INF1.jma > jma
R --no-save -q <<END
  require(gap);require(gap.examples)
  jma <- read.delim("jma", as.is=TRUE)
  primary <- subset(jma,p <= 5e-10 & pJ <= 5e-10)
  secondary <- subset(jma,p > 5e-10 & pJ <= 5e-10)
  print(cbind(nrow(jma),nrow(primary),nrow(secondary)))
  hits <- merge(jma[c("prot","CHR","BP","SNP")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","uniprot")
# total
  cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
  cis.vs.trans <- with(cistrans,data)
  write.table(cis.vs.trans,file="INF1.jma.cis.vs.trans",row.names=FALSE,quote=TRUE)
  cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
  write.table(cis,file="INF1.jma.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
  sink("INF1.jma.out")
  with(cistrans,table)
  sink()
  with(cistrans,total)
  pdf("INF1.jma.pdf")
  circos.cis.vs.trans.plot(hits="jma",inf1,"uniprot")
  dev.off()
# primary
  cistrans1 <- cis.vs.trans.classification(subset(hits,SNP%in%with(primary,SNP)),inf1,"uniprot")
  cis.vs.trans1 <- with(cistrans1,data)
  sink("INF1.jma.1.out")
  with(cistrans1,table)
  sink()
# secondary
  cistrans2 <- cis.vs.trans.classification(subset(hits,SNP%in%with(secondary,SNP)),inf1,"uniprot")
  cis.vs.trans2 <- with(cistrans2,data)
  sink("INF1.jma.2.out")
  with(cistrans2,table)
  sink()
END
rm jma

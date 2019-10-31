# 31-10-2019 JHZ
(
  cat work/*jma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  awk 'NR>1{print $5,$6}' work/INF1.merge | \
  parallel -j1 -C' ' '
    if [ -f work/{1}-{2}.jma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "NR > 1 {print prot, snpid, \$0}" work/{1}-{2}.jma.cojo
    fi'
) > INF1.jma
sed 's/Chr/CHR/g;s/bp/BP/g' INF1.jma > jma
awk 'NR==1 || $10 <= 5e-10' jma > jma.1
awk 'NR==1 || $10 > 5e-10 && $15 <= 5e-10' jma > jma.2
R --no-save -q <<END
  require(gap);require(gap.examples)
  jma <- read.delim("jma", as.is=TRUE)
  primary <- subset(jma,p <= 5e-10)
  secondary <- subset(jma,p > 5e-10 & pJ <= 5e-10)
  print(cbind(nrow(jma),nrow(primary),nrow(secondary)))
  hits <- merge(jma[c("prot","CHR","BP","SNP","p","pJ")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","p","pJ","uniprot")
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
  cistrans1 <- cis.vs.trans.classification(subset(hits,p <= 5e-10),inf1,"uniprot")
  cis.vs.trans1 <- with(cistrans1,data)
  write.table(cis.vs.trans1,file="INF1.jma.1.cis.vs.trans",row.names=FALSE,quote=TRUE)
  sink("INF1.jma.1.out")
  with(cistrans1,table)
  sink()
  with(cistrans1,total)
  pdf("INF1.jma.1.pdf")
  circos.cis.vs.trans.plot(hits="jma.1",inf1,"uniprot")
  dev.off()
# secondary
  cistrans2 <- cis.vs.trans.classification(subset(hits,p > 5e-10 & pJ <= 5e-10),inf1,"uniprot")
  cis.vs.trans2 <- with(cistrans2,data)
  write.table(cis.vs.trans2,file="INF1.jma.2.cis.vs.trans",row.names=FALSE,quote=TRUE)
  sink("INF1.jma.2.out")
  with(cistrans2,table)
  sink()
  with(cistrans2,total)
  pdf("INF1.jma.2.pdf")
  circos.cis.vs.trans.plot(hits="jma.2",inf1,"uniprot")
  dev.off()
END
rm jma jma.1 jma.2
for f in INF1.jma INF1.jma.1 INF1.jma.2
do
  pdftopng -r 300 ${f}.pdf ${f}
  mv ${f}-000001.png ${f}.png
done

R --no-save -q <<END
  library(gap)
  d <- read.table("INF1.jma.cis.vs.trans",as.is=TRUE,header=TRUE)
  pdf("INF1.jma.m2d.pdf")
  mhtplot2d(d)
  dev.off()
  d <- read.table("INF1.jma.1.cis.vs.trans",as.is=TRUE,header=TRUE)
  pdf("INF1.jma.1.m2d.pdf")
  mhtplot2d(d)
  dev.off()
  d <- read.table("INF1.jma.2.cis.vs.trans",as.is=TRUE,header=TRUE)
  pdf("INF1.jma.2.m2d.pdf")
  mhtplot2d(d)
  dev.off()
END

pdftopng -r 300 INF1.jma.m2d.pdf INF1.jma.m2d
mv INF1.jma.m2d-000001.png INF1.jma.m2d.png
pdftopng -r 300 INF1.jma.1.m2d.pdf INF1.jma.1.m2d
mv INF1.jma.1.m2d-000001.png INF1.jma.1.m2d.png
pdftopng -r 300 INF1.jma.2.m2d.pdf INF1.jma.2.m2d
mv INF1.jma.2.m2d-000001.png INF1.jma.2.m2d.png

(
  awk 'NR>1{print $5,$6}' work/INF1.merge | \
  parallel -j1 -C' ' '
    if [ -f work/{1}-{2}.jma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "NR > 1 {print prot \"-\" snpid, NR-1}" work/{1}-{2}.jma.cojo | \
       awk "{a[\$1]=\$2} END {for (i in a) print i, a[i]}"
    fi'
) > work/slct.list

# when slct-INTERVAL.list and slct-ukb.list available
(
  echo locus INTERVAL ukb
  join <(sort -k1,1 slct-INTERVAL.list) <(sort -k1,1 slct-ukb.list) 
) | \
sed 's/ /\t/g' > slct.tsv

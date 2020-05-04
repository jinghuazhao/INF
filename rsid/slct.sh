#!/usr/bin/bash

for j in $(seq 180)
do
  echo ${j}
  sed '1d' ${INF}/work/INF1.merge-rsid | cut -f5,6 | awk -v j=${j} 'NR==j'
  ${INF}/rsid/slct.ini ${j}
done

(
  cat ${INF}/sentinels/*-rsid.jma.cojo | \
  head -1 | \
  awk -v OFS="\t" '{print "prot", "SNPID", $0}'
  awk 'NR>1{print $5,$6}' work/INF1.merge-rsid | \
  parallel -j1 -C' ' '
    if [ -f ${INF}/sentinels/{1}-{2}-rsid.jma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "NR > 1 {print prot, snpid, \$0}" ${INF}/sentinels/{1}-{2}-rsid.jma.cojo
    fi'
) > INF1.jma-rsid
sed 's/Chr/CHR/g;s/bp/BP/g' INF1-rsid.jma > jma
awk 'NR==1 || $10 <= 5e-10' jma > jma.1
awk 'NR==1 || $10 > 5e-10 && $15 <= 5e-10' jma > jma.2
R --no-save -q <<END
  require(gap);require(gap.examples)
  jma <- read.delim("jma", as.is=TRUE)
  primary <- subset(jma,p <= 5e-10)
  secondary <- subset(jma,p > 5e-10 & pJ <= 5e-10)
  print(cbind(nrow(jma),nrow(primary),nrow(secondary)))
  hits <- merge(jma[c("prot","CHR","BP","SNP","b","se","p","bJ","bJ_se","pJ")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","b","se","p","bJ","bJ_se","pJ","uniprot")
# total
  cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
  cis.vs.trans <- with(cistrans,data)
  write.table(cis.vs.trans,file="INF1.jma-rsid.cis.vs.trans",row.names=FALSE,quote=TRUE)
  cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
  write.table(cis,file="INF1.jma-rsid.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
  sink("INF1.jma-rsid.out")
  with(cistrans,table)
  sink()
  with(cistrans,total)
  pdf("INF1.jma-rsid.pdf")
  circos.cis.vs.trans.plot(hits="jma",inf1,"uniprot")
  dev.off()
# primary
  cistrans1 <- cis.vs.trans.classification(subset(hits,p <= 5e-10),inf1,"uniprot")
  cis.vs.trans1 <- with(cistrans1,data)
  write.table(cis.vs.trans1,file="INF1.jma-rsid.1.cis.vs.trans",row.names=FALSE,quote=TRUE)
  sink("INF1.jma-rsid.1.out")
  with(cistrans1,table)
  sink()
  with(cistrans1,total)
  pdf("INF1.jma-rsid.1.pdf")
  circos.cis.vs.trans.plot(hits="jma.1",inf1,"uniprot")
  dev.off()
# secondary
  cistrans2 <- cis.vs.trans.classification(subset(hits,p > 5e-10 & pJ <= 5e-10),inf1,"uniprot")
  cis.vs.trans2 <- with(cistrans2,data)
  write.table(cis.vs.trans2,file="INF1.jma-rsid.2.cis.vs.trans",row.names=FALSE,quote=TRUE)
  sink("INF1.jma-rsid.2.out")
  with(cistrans2,table)
  sink()
  with(cistrans2,total)
  pdf("INF1.jma-rsid.2.pdf")
  circos.cis.vs.trans.plot(hits="jma.2",inf1,"uniprot")
  dev.off()
END
rm jma jma.1 jma.2
for f in INF1.jma-rsid INF1.jma-rsid.1 INF1.jma-rsid.2
do
  pdftopng -r 300 ${f}.pdf ${f}
  mv ${f}-000001.png ${f}.png
done

R --no-save -q <<END
  library(gap)
  d <- read.table("INF1.jma-rsid.cis.vs.trans",as.is=TRUE,header=TRUE)
  d <- within(d,{log10p=-log10p(b/se)})
  pdf("INF1.jma-rsid.m2d.pdf")
  mhtplot2d(d)
  dev.off()
  d <- read.table("INF1.jma-rsid.1.cis.vs.trans",as.is=TRUE,header=TRUE)
  d <- within(d,{log10p=-log10p(b/se)})
  pdf("INF1.jma-rsid.1.m2d.pdf")
  mhtplot2d(d)
  dev.off()
  d <- read.table("INF1.jma-rsid.2.cis.vs.trans",as.is=TRUE,header=TRUE)
  d <- within(d,{log10p=-log10p(b/se)})
  pdf("INF1.jma-rsid.2.m2d.pdf")
  mhtplot2d(d)
  dev.off()
END

pdftopng -r 300 INF1.jma-rsid.m2d.pdf INF1.jma-rsid.m2d
mv INF1.jma-rsid.m2d-000001.png INF1.jma-rsid.m2d.png
pdftopng -r 300 INF1.jma-rsid.1.m2d.pdf INF1.jma-rsid.1.m2d
mv INF1.jma-rsid.1.m2d-000001.png INF1.jma-rsid.1.m2d.png
pdftopng -r 300 INF1.jma-rsid.2.m2d.pdf INF1.jma-rsid.2.m2d
mv INF1.jma-rsid.2.m2d-000001.png INF1.jma-rsid.2.m2d.png

(
  awk 'NR>1{print $5,$6}' ${INF}/work/INF1.merge-rsid | \
  parallel -j1 -C' ' '
    if [ -f sentinels/{1}-{2}-rsid.jma.cojo ]; then
       awk -vprot={1} -v snpid={2} -vOFS="\t" "NR > 1 {print prot \"-\" snpid, NR-1}" ${INF}/sentinels/{1}-{2}-rsid.jma.cojo | \
       awk "{a[\$1]=\$2} END {for (i in a) print i, a[i]}"
    fi'
) > ${INF}/sentinels/slct-rsid.list

# dropped from --cojo-slct
join -v2 <(sed '1d' sentinels/INF1.jma-rsid | cut -f1,2 | tr '\t' '-' | sort) \
         <(sed '1d' work/INF1.merge-rsid | cut -f5,6 | tr '\t' '-' | sort)

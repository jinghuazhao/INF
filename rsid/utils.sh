#!/usr/bin/bash

# count genomic regions for all sentinels
(
  awk -v OFS='\t' '(NR==1){print $1,$2,$3,$5,$6,$8,$9}' work/INF1.merge 
  awk -vd=1e6 -v OFS='\t' '
    (NR>1){
    if($3-$2<=2) {$2=$2-d;$3=$3+d}
    if ($2<0) $2=0
    print $1,$2,$3,$5,$6,$8,$9
  }' work/INF1.merge | \
  sort -k6,6n -k2,2n
) | \
bedtools merge | \
wc -l

awk -vFS="," -vOFS="\t" 'NR>1{print $3,$22,$21}' INF1.jma-rsid.cis.vs.trans | \
sort -k1n -k2n | \
awk '{print "chr" $0}' | \
bedtools merge | \
wc -l

# ST2
R --no-save -q <<END
  cvt<- read.table("work/INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  m <- read.delim("work/INF1.merge",as.is=TRUE)
  m <- merge(m[c("prot","MarkerName","Start","End")],cvt[c("uniprot","prot","SNP","cis.trans")],by.x=c("prot","MarkerName"),by.y=c("prot","SNP"))
  load("ds/latest/fp/INF1.rda")
  e <- with(tbl,prot=="CCL25" & MarkerName=="chr19:49206145_C_G")
  tbl <- subset(tbl,!e)
  m <- merge(tbl,m,by=c("prot","MarkerName"))
  m <- merge(rsid,m,by="MarkerName")
  s <- setdiff(names(m),c("FreqSE","MinFreq","MaxFreq"))
  write.table(m[s],file="work/INF1.METAL",quote=FALSE,row.names=FALSE,sep="\t")
END

# CNVplot
R --no-save -q <<END
  png("work/regions.png",width=10,height=8,units="in",pointsize=8,res=300)
  xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y";
  metal <- read.delim("work/INF1.METAL",as.is=TRUE)[c("Chromosome","Start","End", "prot", "cis.trans")]
  names(metal) <- c("chr","start","end","prot","cis.trans")
  d <- within(metal,{chr<-replace(chr,chr=="X",23); chr<-replace(chr,chr=="Y",24)})
  t <- table(with(metal,chr))
  n <- length(t)
  pos <- vector("numeric")
  for (x in 1:n) pos[x] <- with(subset(d,chr==paste(x)),max(end))
  CM <- cumsum(pos)
  par(xaxt = "n", yaxt = "n")
  xycoords <- xy.coords(c(0,CM), seq(1,5,by=4/n))
  with(xycoords,plot(x, y, type = "n", ann = FALSE, axes = FALSE))
  par(xaxt = "s", yaxt = "s", xpd = TRUE)
  for (x in 1:n) with(subset(d,chr==paste(x)), {
      l <- ifelse(x==1,0,CM[x-1])
      pos[x] <<- ifelse(x == 1, CM[x]/2, (CM[x-1] + CM[x])/2)
      h <- 1+1:5/t[x]
      segments(l+start,h,l+end,h,lwd="3",col=ifelse(cis.trans=="cis","red","blue"))
  })
  axis(1,labels=names(t),at=pos)
  title(main="Flanking regions of sentinels (red=cis, blue=trans)",xlab="Chromosome",ylab="",line=2.5)
  dev.off()
END

R --no-save -q <<END
  cvt <- read.table("work/INF1.merge.out",as.is=TRUE,header=TRUE,nrows=70)
  H <- with(cvt,table(total))
  M <- names(H)
  png(file = "work/signals_by_protein.png",width=6,height=5,units="cm",pointsize=4,res=300)
  barplot(H,names.arg=M,xlab="No. of pQTL regions",ylab="No. of proteins",ylim=c(0,25),col="darkgrey",border="black")
  dev.off()
END

# twas_fusion coloc results
function get_fusion()
{
  export tissue=Whole_Blood
  export tag=${tissue}-coloc.dat
  cd ${INF}/work/coloc
  (
    cat *${tag} | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tag} | awk -vtag=-${tag} '{sub(tag,"")};1' | \
    parallel --env tag -C' ' 'grep -v PP4 {}-${tag} | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fusion
}

# fastenloc coloc results
function get_fastenloc()
{
  export tissue=Whole_Blood
  cd ${INF}/work/fastenloc
  (
    cat *${tissue}.sig.out | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tissue}.sig.out | awk -vtag=-Whole_Blood.sig.out '{sub(tag,"")};1' | \
    parallel -C' ' 'grep -v RCP {}-${tissue}.sig.out | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fastenloc-sig
  (
    cat *${tissue}.snp.out | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tissue}.snp.out | awk -vtag=-${tissue}.snp.out '{sub(tag,"")};1' | \
    parallel -C' ' 'grep -v SCP {}-${tissue}.snp.out | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fastenloc-snp
}

# eQTL overlap
awk '!/BDNF/ && NR > 1 {
  if ($3=="\"Q8NF90\"") $7="\"FGF5\""; else if ($3=="\"Q8WWJ7\"") $7="\"CD6\"";print}
' FS='\t' OFS='\t' doc/olink.inf.panel.annot.tsv | \
cut -f2,3,7 | \
cut -f3 | \
sed 's/"//g' > work/inf1.gene

sed '1d' work/INF1.merge.cis.vs.trans | \
cut -d, -f10 | \
uniq > work/INF1.gene

# phenoscanner -t T -c eQTL -x EUR -p 5e-8 -r 0.8 -i work/INF1.gene -o INF1
# additionally from pqtlxQTL.R

function ieu()
{
  md5sum * > INF1.md5sum
  md5sum --check INF1.md5sum
  ls  *gz | parallel -C' ' 'echo {}; gunzip -c {} | wc -l' | paste - - | sort -k1,1 > INF1.size
}

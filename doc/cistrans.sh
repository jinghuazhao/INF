#!/bin/bash
. /etc/profile.d/modules.sh

R --no-save -q <<END

# preliminary annotation
inf1 <- read.delim("/scratch/jhz22/INF/doc/olink.inf.panel.annot.tsv", as.is=TRUE)
inf1[with(inf1, uniprot=="Q8NF90"),"hgnc_symbol"] <- "FGF5"
inf1[with(inf1, uniprot=="Q8WWJ7"),"hgnc_symbol"] <- "CD6"

prot <- read.table("/scratch/jhz22/INF/inf1.list",col.names=c("prot","uniprot"),as.is=TRUE,sep="\t")

p <- merge(inf1,prot,by="uniprot")[c("chromosome_name","start_position","end_position","hgnc_symbol","prot")]
p <- within(p,{chromosome_name=paste0("chr",chromosome_name)})
names(p) <- c("#chrom","start","end","gene","prot")
write.table(p,file="inf1.bed",quote=FALSE,row.names=FALSE,sep="\t")

jma <- read.table("INTERVAL.jma.dat",as.is=TRUE,header=TRUE)
r <- merge(jma[c("Chr","bp","SNP","prot")],p[c("prot","gene")],by="prot")
write.table(r,file="cistrans.tmp",quote=FALSE,row.names=FALSE,sep="\t")
END

awk -vOFS="\t" '{
  if(NR==1) print "#chrom","start","end","snp","prot","gene";
  else print "chr" $2, $3-1, $3, $4, $1, $5
}' cistrans.tmp > cistrans.in
awk -vOFS="\t" -vM=1000000 '{
  if(NR==1) print;
  else {
    start=$2-M;
    if(start<0) start=0;
    end=$3+M;
    print $1,start,end,$4,$5
  }
}' inf1.bed > inf1.tmp
module load gcc/4.8.1
(
  echo -e "chr\tstart\tend\tSNP\tprot\tgene\tstatus"
  bedtools intersect -a cistrans.in -b inf1.tmp -u > cistrans.tmp
  bedtools intersect -a cistrans.in -b inf1.tmp -loj -u > cistrans.tmp
  bedtools intersect -a cistrans.in -b inf1.tmp -loj | \
  awk '($6==$10)' | \
  cut -f4 > cistrans.snp
  grep -w -f cistrans.snp cistrans.tmp | \
  awk -vOFS="\t" '{print $0, "cis"}'
  grep -v -w -f cistrans.snp cistrans.tmp | \
  awk -vOFS="\t" '{print $0, "trans"}'
  bedtools intersect -a cistrans.in -b inf1.tmp -v | \
  awk -vOFS="\t" '{print $0,"trans"}'
) > cistrans.tsv

(
R --no-save -q <<END
  cistrans <- read.delim("cistrans.tsv",as.is=TRUE)
  with(cistrans, table(gene,status))
END
) > cistrans.table

#!/bin/bash
. /etc/profile.d/modules.sh

export rt=/scratch/jhz22/INF
export prot_list=$rt/doc/olink.prot.list.txt
export prot_annotation=$rt/doc/olink.inf.panel.annot.tsv
export prot_jma_cojo=INTERVAL.jma.dat
(
  grep inf1 ${prot_list} | \
  sed 's/inf1_//g;s/___/\t/g'
) | \
sort -k1,1 > inf1.list
(
  echo -e "chrom\tstart\tend\tgene\tprot"
  sort -k2,2 inf1.list > inf1.tmp
  cut -f2,3,7-10 ${prot_annotation}  | \
  awk -vFS="\t" -vOFS="\t" '(NR>1){
      gsub(/\"/,"",$0)
      if($2=="Q8NF90") $3="FGF5"
      if($2=="Q8WWJ7") $3="CD6"
      print
  }' | \
  sort -t$'\t' -k2,2 | \
  join -t$'\t' -j2 inf1.tmp - | \
  awk -vFS="\t" -vOFS="\t" '{print $5,$6,$7,$4,$2}' | \
  sort -k1,1n -k2,2n | \
  awk -vFS="\t" -vOFS="\t" '{$1="chr" $1;print}'
) > inf1.bed
(
  echo -e "prot\tChr\tbp\tSNP\tgene"
  awk -vFS="\t" -vOFS="\t" 'NR>1 {print $4,$5}' inf1.bed | \
  sort -k2,2 > inf1.tmp
  awk -vOFS="\t" 'NR>1 {print $1,$2,$3,$4}' ${prot_jma_cojo} | \
  sort -k1,1 | \
  join -t$'\t' -11 -22 - inf1.tmp | \
  awk -vFS="\t" -vOFS="\t" '{print $1,$2,$4,$3,$5}'
) > cistrans.tmp

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
  bedtools intersect -a cistrans.in -b inf1.tmp -loj | \
  awk '($6==$10)' | \
  cut -f4 > cistrans.cis
  grep -w -f cistrans.cis cistrans.tmp | \
  awk -vOFS="\t" '{print $0, "cis"}'
  grep -v -w -f cistrans.cis cistrans.tmp | \
  awk -vOFS="\t" '{print $0, "trans"}'
  bedtools intersect -a cistrans.in -b inf1.tmp -v | \
  awk -vOFS="\t" '{print $0,"trans"}'
) > cistrans.tsv

R --no-save -q <<END
  sink("cistrans.table")
  cistrans <- read.delim("cistrans.tsv",as.is=TRUE)
  with(cistrans, table(gene,status))
  sink()

# for R, data preparation is through the following section,
# while the bedtools section can be furnished with Bioconductor/GenomicRanges

  preproc <- function()
  {
    inf1 <- read.delim("olink.inf.panel.annot.tsv", as.is=TRUE)
    inf1[with(inf1, uniprot=="Q8NF90"),"hgnc_symbol"] <- "FGF5"
    inf1[with(inf1, uniprot=="Q8WWJ7"),"hgnc_symbol"] <- "CD6"

    prot <- read.table("inf1.list",col.names=c("prot","uniprot"),as.is=TRUE,sep="\t")

    p <- merge(inf1,prot,by="uniprot")[c("chromosome_name","start_position","end_position","hgnc_symbol","prot")]
    p <- within(p,{chromosome_name=paste0("chr",chromosome_name)})
    names(p) <- c("#chrom","start","end","gene","prot")
    write.table(p,file="inf1.bed",quote=FALSE,row.names=FALSE,sep="\t")

    jma <- read.table("INTERVAL.jma.dat",as.is=TRUE,header=TRUE)
    r <- merge(jma[c("Chr","bp","SNP","prot")],p[c("prot","gene")],by="prot")
    write.table(r,file="cistrans.tmp",quote=FALSE,row.names=FALSE,sep="\t")
  }
END

cat cistrans.table

# head -1 $rt/doc/olink.inf.panel.annot.tsv |  awk '{gsub(/\t/, "\n",$0)};1'|  awk '{print "#", NR, $1}'
# 1 "target"
# 2 "target.short"
# 3 "uniprot"
# 4 "panel"
# 5 "prot.on.multiple.panel"
# 6 "panels.with.prot"
# 7 "hgnc_symbol"
# 8 "chromosome_name"
# 9 "start_position"
# 10 "end_position"
# 11 "olink.id"
# 12 "alternate.uniprot"

# head -1 INTERVAL.jma.dat | awk '{gsub(/ /, "\n",$0);print}' | awk '{print "#", NR, $1}'
# 1 prot
# 2 Chr
# 3 SNP
# 4 bp
# 5 refA
# 6 freq
# 7 b
# 8 se
# 9 p
# 10 n
# 11 freq_geno
# 12 bJ
# 13 bJ_se
# 14 pJ
# 15 LD_r

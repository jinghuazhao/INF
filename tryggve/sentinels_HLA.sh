# 6-6-2019 JHZ

module load bedtools/2.27.1

function pgz()
# extract all significant SNPs to .p.gz
{
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j3 -C' ' '
  (
  # zcat METAL/{}-1.tbl.gz | head -1
    zcat METAL/{}-1.tbl.gz | awk "NR>1 && length(\$4)==1 && length(\$5)==1 && \$12<5e-10" | sort -k1,1n -k2,2n
  ) | gzip -f > work/{}.p.gz'
}

function nold_HLA()
# removing those in high LD regions'
{
  awk '($4!=8)' tryggve/high-LD-regions-hg19.bed > work/high-LD-regions-HLA-hg19.bed
  for p in $(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g')
  do
    (
      zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      zcat work/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}'
    ) | bedtools subtract -header -a - -b work/high-LD-regions-HLA-hg19.bed > work/${p}.p
    (
      head -1 work/${p}.p
      export lines=$(wc -l work/${p}.p | cut -d' ' -f1)
      if [ $lines -gt 1 ]; then
        awk 'NR>1' work/${p}.p
        awk '$1 == "chr6" && $3 >= 25392021 && $3 < 33392022' work/${p}.p | \
        sort -k13,13g | \
        awk 'NR==1'
      fi
    ) > work/${p}_HLA.p
    rm work/${p}.p
  done
}

function sentinels()
# find sentinels
{
  for prot in $(ls work/*_HLA.p | sed 's|work/||g;s|_HLA.p||g')
  do 
    export prot=${prot}
    export tag=_HLA
    echo ${prot}${tag}
    R --no-save -q < tryggve/sentinels.R > work/${prot}_HLA.o
  done
  cd work
  (
    awk -vOFS="," 'BEGIN{print "prot","CHR","BP","SNP","l","u","d","log10p","Groupid", "Type"}'
    awk -vFS="," -vOFS="," '!/option/{
        SNPID=$2
        split(SNPID,a,":")
        split(a[2],b,"_")
        gsub(/chr/,"",a[1])
        $1=$1 "," a[1] "," b[1]
        print
    }' *_HLA.o
  ) | sed 's/,/ /g' > INF1_HLA.sentinels
  cd -
}

function cvt()
# cis.vs.classification, requiring R-3.5.0 at TRYGGVE
{
  cd work
  R --no-save -q <<\ \ END
    require(gap)
    clumped <- read.table("INF1_HLA.sentinels",as.is=TRUE,header=TRUE)
    hits <- merge(clumped[c("CHR","BP","SNP","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
    cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
    cis.vs.trans <- with(cistrans,data)
    write.table(cis.vs.trans,file="INF1_HLA.sentinels.cis.vs.trans",row.names=FALSE,quote=TRUE)
    cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
    write.table(cis,file="INF1_HLA.sentinels.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
    sink("INF1_HLA.sentinels.out")
    with(cistrans,table)
    sink()
    with(cistrans,total)
    pdf("INF1_HLA.sentinels.circlize.pdf")
    circos.cis.vs.trans.plot(hits="INF1_HLA.sentinels",inf1,"uniprot")
    dev.off()
  END
  cd -
}

$1

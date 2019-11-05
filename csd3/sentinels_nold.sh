# 5-11-2019 JHZ

export tag=_nold

function pgz()
# 1. extract all significant SNPs
{
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j3 -C' ' '
  (
  # zcat METAL/{}-1.tbl.gz | head -1
    zcat METAL/{}-1.tbl.gz | awk "
    function abs(x)
    {
      if (x<0) return -x;
      else return x;
    }
    NR>1 && length(\$4)==1 && length(\$5)==1 && abs(\$10/\$11)>=6.219105" | sort -k1,1n -k2,2n
  ) | gzip -f > sentinels/{}.p.gz'
}

function _HLA()
# 2. handling HLA
{
  for p in $(ls METAL/*-1.tbl.gz | sed 's|METAL/||g;s/-1.tbl.gz//g')
  do
    (
      zcat METAL/${p}-1.tbl.gz | head -1 | awk -vOFS="\t" '{$1="Chrom";$2="Start" "\t" "End";print}'
      zcat sentinels/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
      awk '!($1 == "chr6" && $3 >= 25392021 && $3 < 33392022)'
      zcat sentinels/${p}.p.gz | \
      awk -vOFS="\t" '{$1="chr" $1; start=$2-1;$2=start "\t" $2;print}' | \
      awk '$1 == "chr6" && $3 >= 25392021 && $3 < 33392022' | \
      sort -k13,13g | \
      awk 'NR==1'
    ) > sentinels/${p}${tag}.p
    export lines=$(wc -l sentinels/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm sentinels/${p}${tag}.p
    fi
  done
}

function sentinels()
# 3. find sentinels
{
  for prot in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g')
  do 
    export prot=${prot}
    echo ${prot}${tag}
    R --no-save -q < tryggve/sentinels.R > sentinels/${prot}${tag}.o
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
    }' *${tag}.o
  ) | sed 's/,/ /g' > INF1${tag}.sentinels
  awk '$1 != "total" && NR > 1' INF1${tag}.sentinels.out | wc -l
  awk '$1 != "total" && NR > 1 && $2 > 0 && $3 == 0' INF1${tag}.sentinels.out | wc -l
  awk '$1 != "total" && NR > 1 && $2 == 0 && $3 > 0' INF1${tag}.sentinels.out | wc -l
  awk '$1 != "total" && NR > 1 && $2 > 0 && $3 > 0' INF1${tag}.sentinels.out | wc -l
  cd -
}

function cvt()
# 4. cis.vs.classification, requiring R-3.5.0 at TRYGGVE
{
  cd work
  R --no-save -q <<\ \ END
    require(gap)
    tag <- Sys.getenv("tag")
    rt <- paste0("INF1",tag)
    clumped <- read.table(paste0(rt,".sentinels"),as.is=TRUE,header=TRUE)
    hits <- merge(clumped[c("CHR","BP","SNP","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
    cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
    cis.vs.trans <- with(cistrans,data)
    write.table(cis.vs.trans,file=paste0(rt,".sentinels.cis.vs.trans"),row.names=FALSE,quote=TRUE)
    cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
    write.table(cis,file=paste0(rt,".sentinels.cis"),col.names=FALSE,row.names=FALSE,quote=FALSE)
    sink(paste0(rt,".sentinels.out"))
    with(cistrans,table)
    sink()
    with(cistrans,total)
    pdf(paste0(rt,".sentinels.circlize.pdf"))
    circos.cis.vs.trans.plot(hits=paste0(rt,".sentinels"),inf1,"uniprot")
    dev.off()
  END
  cd -
}

for cmd in pgz _HLA sentinels cvt; do $cmd; done

# 18-6-2019 JHZ

module load bedtools/2.27.1
export tag=_nold

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
    ) > work/${p}${tag}.p
    export lines=$(wc -l work/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm work/${p}${tag}.p
    fi
  done
}

function sentinels()
# 3. find sentinels
{
  for prot in $(ls work/*${tag}.p | sed 's|work/||g;s|'"$tag"'.p||g')
  do 
    export prot=${prot}
    echo ${prot}${tag}
    R --no-save -q < tryggve/sentinels.R > work/${prot}${tag}.o
  done
  cd work
  (
    awk -vOFS="," 'BEGIN{print "prot","CHR","BP","SNP","l","u","d","logp","Groupid", "Type"}'
    awk -vFS="," -vOFS="," '!/option/{
        SNPID=$2
        split(SNPID,a,":")
        split(a[2],b,"_")
        gsub(/chr/,"",a[1])
        $1=$1 "," a[1] "," b[1]
        print
    }' *${tag}.o
  ) | sed 's/,/ /g' > INF1${tag}.sentinels
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
    hits <- merge(clumped[c("CHR","BP","SNP","prot","logp")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","logp","uniprot")
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

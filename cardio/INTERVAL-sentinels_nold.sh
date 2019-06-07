# 7-6-2019 JHZ

export tag=_nold

function pgz()
# 1. extract all significant SNPs to .p.gz
{
  ls sumstats/INTERVAL/INTERVAL*.gz | \
  xargs -l basename | \
  sed 's|INTERVAL.||g;s/.gz//g' | \
  parallel -j3 -C' ' '
  (
    zcat sumstats/INTERVAL/INTERVAL.{}.gz | awk "NR>1 && length(\$6)==1 && length(\$7)==1 && \$11<5e-10" | sort -k2,2n -k3,3n
  ) | gzip -f > INTERVAL/{}.p.gz'
}

function _HLA()
# 2. handling HLA
{
  for p in $(ls sumstats/INTERVAL/INTERVAL*.gz | \
  xargs -l basename | \
  sed 's|INTERVAL.||g;s/.gz//g')
  do
    (
      zcat sumstats/INTERVAL/INTERVAL.${p}.gz | head -1 | \
      awk -vOFS="\t" '{$1="MarkerName";$2="Chrom";$3="Start" "\t" "End";print}'
      zcat INTERVAL/${p}.p.gz | \
      awk -vOFS="\t" '{$2="chr" $2; start=$3-1;$3=start "\t" $3;print}' | \
      awk '!($2 == "chr6" && $4 >= 25392021 && $4 < 33392022)'
      zcat INTERVAL/${p}.p.gz | \
      awk -vOFS="\t" '{$2="chr" $2; start=$3-1;$3=start "\t" $3;print}' | \
      awk '$2 == "chr6" && $4 >= 25392021 && $4 < 33392022' | \
      sort -k12,12g | \
      awk 'NR==1'
    ) > INTERVAL/${p}${tag}.p
    export lines=$(wc -l INTERVAL/${p}${tag}.p | cut -d' ' -f1)
    if [ $lines -eq 1 ]; then
      echo removing ${p}${tag} with $lines lines
      rm INTERVAL/${p}${tag}.p
    fi
  done
}

function sentinels()
# 3. find sentinels
{
  for prot in $(ls INTERVAL/*${tag}.p | sed 's|INTERVAL/||g;s|'"$tag"'.p||g')
  do 
    export prot=${prot}
    echo ${prot}${tag}
    R --no-save -q < cardio/sentinels.R > INTERVAL/${prot}${tag}.o
  done
  cd INTERVAL
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
  ) | sed 's/,/ /g' > INTERVAL${tag}.sentinels
  cd -
}

function cvt()
# 4. cis.vs.classification
{
  cd work
  R --no-save -q <<\ \ END
    require(gap)
    tag <- Sys.getenv("tag")
    rt <- paste0("INTERVAL",tag)
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

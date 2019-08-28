# 28-8-2019 JHZ

source tryggve/analysis.ini

gunzip -c sumstats/INTERVAL/INTERVAL.4E.BP1.gz | head -1 > work/sumstats.hdr
gunzip -c METAL/4E.BP1-1.tbl.gz | head -1 > work/METAL.hdr

function qml()
{
  echo "--> Q-Q/Manhattan/LocusZoom plots"
  export rt=$HOME/INF
  (
    echo -e "chrom\tstart\tend\tgene\tprot"
    sort -k2,2 $rt/inf1.list > inf1.tmp
    cut -f2,3,7-10 $rt/doc/olink.inf.panel.annot.tsv  | \
    awk -vFS="\t" -vOFS="\t" '(NR>1){
        gsub(/\"/,"",$0)
        if($2=="Q8NF90") $3="FGF5"
        if($2=="Q8WWJ7") $3="CD6"
        print
    }' | \
    sort -t$'\t' -k2,2 | \
    join -t$'\t' -j2 inf1.tmp - | \
    awk -vFS="\t" -vOFS="\t" '{print $5,$6,$7,$4,$2}' | \
    sort -k1,1n -k2,2n
  ) > st.bed
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j4 --env rt -C' ' 'export protein={}; R --no-save -q < $rt/tryggve/qqman.R'
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j1 -C' ' '
  (
     echo -e "MarkerName\tP-value\tWeight"
     grep -w {} st.bed | \
     awk -vOFS="\t" -vM=1000000 "{start=\$2-M;if(start<0) start=0;end=\$3+M;\$2=start;\$3=end};1" > st.tmp
     read chrom start end gene prot < st.tmp
     gunzip -c METAL/{}-1.tbl.gz | \
     awk -vOFS="\t" -vchr=$chrom -vstart=$start -vend=$end \
         "(\$1 == chr && \$2 >= start && \$2 <= end){split(\$3,a,\"_\");print a[1],10^\$12,\$18}" | \
     sort -k1,1 | \
     join -12 -21 work/snp_pos - | \
     awk -vOFS="\t" "{print \$2, \$3, \$4}"
  )  > METAL/{}.lz'
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j1 -C' ' '
     grep -w {} st.bed | \
     awk -vOFS="\t" -vM=1000000 "{start=\$2-M;if(start<0) start=0;end=\$3+M;\$2=start;\$3=end};1" > st.tmp
     read chrom start end gene prot < st.tmp
     cd METAL
     rm -f ld_cache.db
     locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal {}.lz \
               --plotonly --chr $chrom --start $start --end $end --no-date --rundir .
     mv chr${chrom}_${start}-${end}.pdf {}.lz.pdf
     pdftopng -r 300 {}.lz.pdf {}
     mv {}-000001.png {}.lz-1.png
     mv {}-000002.png {}.lz-2.png
     cd -
  '
  export PATH=/data/jinhua/ImageMagick-7.0.8-22/bin:%PATH
  convert OPG.lz-1.png -resize 130% OPG.lz-3.png
  convert \( OPG.qq.png -append OPG.manhattan.png -append OPG.lz-3.png -append \) +append OPG-qml.png
#  convert OPG.lz-1.png -resize 130% OPG.lz-3.png
#  pdftopng -r 300 -f 310 -l 310 INF1.fp.pdf 310
#  convert 310-000310.png -resize 150% 310.png
#  convert \( OPG.qq.png -append OPG.manhattan.png -append OPG.lz-3.png -append 310.png -append \) +append OPG-qml.png
}

function clumping()
{
  echo "--> clumping"
  export rt=$HOME/INF/METAL
  ls METAL/*tbl.gz | \
  sed 's/-1.tbl.gz//g' | \
  xargs -l basename | \
  parallel -j4 --env rt -C' ' '
  if [ -f $rt/{}.clumped ]; then rm $rt/{}.clumped; fi
  plink --bfile EUR \
        --clump $rt/{}-1.tbl.gz \
        --clump-snp-field MarkerName \
        --clump-field P-value \
        --clump-kb 1000 \
        --clump-p1 5e-10 \
        --clump-p2 0.01 \
        --clump-r2 0.1 \
        --mac 50 \
        --out $rt/{}
  '
  (
    grep CHR $rt/*.clumped | \
    head -1
    ls METAL/*-1.tbl.gz | \
    sed 's|METAL/||g;s/-1.tbl.gz//g' | \
    parallel -j1 --env rt -C' ' 'grep -H -v CHR $rt/{}.clumped | sed "s/.clumped://g"'
  ) | \
  sed 's|'"$rt"'/||g;s/.clumped://g' | \
  awk '(NF>1){$3="";print}' | \
  awk '{$1=$1;if(NR==1)$1="prot";print}' > INF1.clumped
  R --no-save -q <<\ \ END
    require(gap); require(gap.examples)
    clumped <- read.table("INF1.clumped",as.is=TRUE,header=TRUE)
    hits <- merge(clumped[c("CHR","BP","SNP","prot")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","uniprot")
    cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
    cis.vs.trans <- with(cistrans,data)
    write.table(cis.vs.trans,file="INF1.clumped.cis.vs.trans",row.names=FALSE,quote=TRUE)
    cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
    write.table(cis,file="INF1.clumped.cis",col.names=FALSE,row.names=FALSE,quote=FALSE)
    sink("INF1.clumped.out")
    with(cistrans,table)
    sink()
    with(cistrans,total)
    pdf("INF1.circlize.pdf")
    circos.cis.vs.trans.plot(hits="INF1.clumped",inf1,"uniprot")
    dev.off()
  END
}

function ma()
{
  echo "--> .ma"
  export rt=$HOME/INF/METAL
  ls $rt/*.tbl.gz | \
  sed 's/-1.tbl.gz\*//g' | \
  xargs -l basename | \
  parallel -j4 --env rt -C' ' '
  (
    echo SNP A1 A2 freq b se p N;
    gunzip -c $rt/{}-1.tbl.gz | \
    awk "(NR>1 && \$18>50) {print \$3, toupper(\$4), toupper(\$5), \$6, \$10, \$11, 10^\$12, \$18}"
  ) > work/{}.ma
  '
}
sed 's/\t/\n/g' work/METAL.hdr| awk '{print "# " NR " " $1}'
# 1 Chromosome
# 2 Position
# 3 MarkerName
# 4 Allele1
# 5 Allele2
# 6 Freq1
# 7 FreqSE
# 8 MinFreq
# 9 MaxFreq
# 10 Effect
# 11 StdErr
# 12 log(P)
# 13 Direction
# 14 HetISq
# 15 HetChiSq
# 16 HetDf
# 17 logHetP
# 18 N

function cojo()
{
  echo "--> COJO analysis"
  export rt=$HOME/INF/METAL
  ls METAL/*.tbl.gz | \
  sed 's/-1.tbl.gz//g' | \
  xargs -l basename | \
  parallel -j2 --env rt -C' ' '
  if [ -f $rt/{}.jma.cojo ]; then rm $rt/{}.jma.cojo $rt/{}.ldr.cojo; fi; \
  gcta64 --bfile EUR --cojo-file $rt/{}.ma --cojo-slct --cojo-p 5e-10 --cojo-collinear 0.1 --cojo-wind 1000 \
         --maf 0.01 --thread-num 2 --out $rt/{}
  '
  (
    grep SNP $rt/*.jma.cojo | \
    head -1
    grep -v -w SNP $rt/*.jma.cojo
  ) | \
  sed 's|'"$rt"'/||g;s/.jma.cojo:/\t/g' | \
  awk -vOFS="\t" '{if(NR==1) $1="prot";print}' > INF1.jma
  sed 's/Chr/CHR/g;s/bp/BP/g' INF1.jma > jma
  R --no-save -q <<\ \ END
    require(gap);require(gap.examples)
    jma <- read.delim("jma")
    ni <- dim(jma)[1]
    primary <- dim(subset(jma,p <= 5e-10))[1]
    secondary <- dim(subset(jma,p > 5e-10 & pJ <= 5e-10))[1]
    print(cbind(ni,primary,secondary))
    hits <- merge(jma[c("prot","CHR","BP","SNP")],inf1[c("prot","uniprot")],by="prot")
    names(hits) <- c("prot","Chr","bp","SNP","uniprot")
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
  END
  rm jma
}
#1 MarkerName
#2 Allele1
#3 Allele2
#4 Freq1
#5 FreqSE
#6 MinFreq
#7 MaxFreq
#8 Effect
#9 StdErr
#10 P-value
#11 Direction
#12 CHR
#13 POS
#14 WEIGHT

function fp()
{
  (
    cat work/METAL.hdr
  # 1. PLINK --clumping. The default.
  # 3. GCTA --cojo.  Change Chr, Pos to CHR, BP as in 1.
  # 2. R/gap/sentinels output. use  SNPID or $4 below
    awk 'NR>1 {print $1,$4}' work/INF1_nold.sentinels | \
    parallel -j4 -C' ' 'zgrep -w -H {2} METAL/{1}-1.tbl.gz'
  ) | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' > INF1.tbl
  cut -f3 INF1.tbl | \
  awk 'NR>1' | \
  sort -k1,1 | \
  uniq | \
  join -a2 -e "NA" work/INTERVAL.rsid - -o2.1,1.2> INF1.rsid
  (
    cat work/sumstats.hdr
    awk 'NR>1' INF1.tbl | \
    cut -f1,3,13 | \
    awk '{split($1,a,":");print a[1],$2,$3}' | \
    parallel -j4 -C' ' '
      export direction=$(zgrep -w {2} METAL/{1}-1.tbl.gz | cut -f13)
      let j=1
      for i in $(grep "Input File" METAL/{1}-1.tbl.info | cut -d" " -f7)
      do
         export n=$(awk -vj=$j "BEGIN{split(ENVIRON[\"direction\"],a,\"\");print a[j]}")
         if [ "$n" != "?" ]; then zgrep -H -w {2} $i; fi
         let j=$j+1
      done
  '
  ) | \
  sed 's|/data/jinhua/INF/sumstats||g;s/.gz//g' > INF1.all
  if [ -f INF1.fp.log ]; then rm INF1.fp.log; fi
  script INF1.fp.log
  R -q --no-save <<\ \ END
    require(gap)
    t <- read.delim("INF1.tbl",as.is=TRUE)
    tbl <- within(t, {
      prot <- sapply(strsplit(Chromosome,":"),"[",1)
      Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    })
    a <- read.table("INF1.all",as.is=TRUE, header=TRUE)
    all <- within(a, {
      dir.study.prot <- sapply(strsplit(SNPID,":"),"[",1)
      p1 <- sapply(strsplit(SNPID,":"),"[",2)
      p2 <- sapply(strsplit(SNPID,":"),"[",3)
      MarkerName <- paste(p1,p2,sep=":")
      study <- sapply(strsplit(dir.study.prot,"/"),"[",2)
      study.prot <- sapply(strsplit(dir.study.prot,"/"),"[",3)
      substudy <- sapply(strsplit(study.prot,"[.]"),"[",1)
      pos <- unlist(lapply(gregexpr("[.]",study.prot),"[",1))
      prot <- substring(study.prot,pos+1)
    })
    droplist <- c("SNPID","dir.study.prot","p1","p2","pos","study.prot","substudy")
    all <- all[setdiff(names(all),droplist)]
    rsid <- read.table("INF1.rsid",as.is=TRUE,col.names=c("MarkerName","rsid"))
    save(tbl,all,rsid,file="INF1.rda",version=2)
    METAL_forestplot(tbl,all,rsid,"INF1.fp.pdf",width=8.75,height=5)
  END
  exit
  grep I2 INF1.fp.log | awk '{gsub(/prot|=|MarkerName|Q|df|p|lower.I2|upper.I2|I2/,"");print}' > INF1.Q
}

function aild()
{
  echo "--> approximately independent LD blocks"
  awk 'NR>1{gsub(/chr/,"",$1);print}' tryggve/EURLD.bed > rlist-EURLD
  export rt=$HOME/INF/AILD
  ls $rt/METAL/*tbl.gz | \
  sed 's/-1.tbl.gz//g' | \
  xargs -l basename | \
  parallel -j5 --env rt -C' ' '
  plink --bfile EUR \
        --clump $rt/METAL/{}-1.tbl.gz --clump-range rlist-EURLD --clump-range-border 250 \
        --clump-snp-field MarkerName \
        --clump-field P-value \
        --clump-p1 5e-10 --clump-p2 0.01 --clump-r2 0.1 \
        --mac 50 \
        --out $rt/{}'
  (
    grep CHR $rt/*.clumped.ranges | \
    head -1
    ls METAL/*-1.tbl.gz | \
    sed 's|METAL/||g;s/-1.tbl.gz//g' | \
    parallel -j1 --env rt -C' ' 'grep -H -v CHR $rt/{}.clumped.ranges | sed "s/.clumped.ranges://g"'
  ) | \
  sed 's|'"$rt"'/||g;s/.clumped://g' | \
  awk '(NF>1){$3="";print}' | \
  awk '{$1=$1;if(NR==1)$1="prot";print}' > INF1.ranges
  awk '(NR>1){
    chr=$1;
    gsub(/chr/,"",chr);
    flanking=($3-$2)/2/1000
    centre=$2+flanking
    print sprintf("%d %d %d %s", chr, centre, flanking, $4)
  }' tryggve/EURLD.bed > rlist=EURLD.region

  for prot in $(ls $rt/METAL/*tbl.gz | sed 's/-1.tbl.gz//g' | xargs -l basename)
  do
    export p=$prot
    # bruteforceclumpingbyregion
    cat rlist=EURLD.region | \
    parallel -j8 --env p --env rt -C' ' '
     gcta64 --bfile EUR --cojo-file $rt/METAL/$p.ma --cojo-slct --cojo-p 5e-10 --maf 0.0001 \
            --extract-region-bp {1} {2} {3} --thread-num 3 --out $rt/LDBLOCK/$p-{4}'
     (
       cat $rt/LDBLOCK/${p}*.jma.cojo | head -1
       awk "NR>1" $rt/LDBLOCK/${p}*.jma.cojo
     ) > $rt/${p}.jma
  done
}

function cs_jma()
{
R -q --no-save <<END
  require(gtx)
  max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
  t <- read.delim("INF1.jma",as.is=TRUE)
  tbl <- within(t, {
    BF.normal <- max1(abf.Wakefield(b, se, 5e-10))
    BF.numeric <- max1(abf.normal(b, se, 5e-10))
    BF.cs <- credset(BF.numeric)
  })
  tbl[c("prot","Chr","p","BF.numeric","BF.cs")]
END
}

function cs_tbl()
{
R -q --no-save <<END
  require(gtx)
  max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
  t <- read.delim("INF1.tbl",as.is=TRUE)
  tbl <- within(t, {
    prot <- sapply(strsplit(Chromosome,":"),"[",1)
    Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    BF.normal <- max1(abf.Wakefield(Effect, StdErr, 0.05))
    BF.numeric <- max1(abf.normal(Effect, StdErr, 0.05))
    BF.cs <- credset(BF.numeric)
  })
  tbl[c("prot","Chromosome","P.value","BF.numeric","BF.cs")]
END
}

function pve()
{
R -q --no-save -q <<END
  t <- read.delim("INF1.tbl",as.is=TRUE)
  tbl <- within(t, {
    prot <- sapply(strsplit(Chromosome,":"),"[",1)
    Chromosome <- sapply(strsplit(Chromosome,":"),"[",2)
    chi2n <- (Effect/StdErr)^2/N
  })
  s <- with(tbl, aggregate(chi2n,list(prot),sum))
  names(s) <- c("prot", "h2")
  sd <- with(tbl, aggregate(chi2n,list(prot),sd))
  names(sd) <- c("p1","sd")
  m <- with(tbl, aggregate(chi2n,list(prot),length))
  names(m) <- c("p2","m")
  h2 <- cbind(s,sd,m)
  ord <- with(h2, order(h2))
  sink("h2.dat")
  print(h2[ord, c("prot","h2","sd","m")], row.names=FALSE)
  sink()
  png("pve.png", res=300, units="in", width=6, height=4)
  np <- nrow(h2)
  with(h2[ord,], {
    plot(h2, cex=0.4, pch=16, xaxt="n", xlab="protein", ylab=expression(h^2))
    xtick <- seq(1, np, by=1)
    axis(side=1, at=xtick, labels = FALSE, lwd.tick=0.2)
    text(x=xtick, par("usr")[3],labels = prot, srt = 75, pos = 1, xpd = TRUE, cex=0.3)
  })
  title("(b)")
  dev.off()
  write.csv(tbl,file="INF1.csv",quote=FALSE,row.names=FALSE)
END
}

function mpfr()
{
  echo "--> METAL results containing P-value=0"
  awk '($5==0)' INF1.clumped | \
  cut -d' ' -f1 | 
  uniq > INF1.z
  cat INF1.z | \
  parall -j2 --env rt '
  (
    export port={}
    gunzip -c $rt/{}-1.tbl.gz | \
    awk -vOFS="\t" "NR==1||\$12!=0"
    gunzip -c $rt/{}-1.tbl.gz | \
    awk -vOFS="\t" "(NR==1||\$12==0)" > {}.z
    R --no-save -q <<\ \ \ \ END
      prot <- Sys.getenv("prot")
      metal <- read.delim(paste0(prot,".z"),as.is=TRUE)
      library(Rmpfr)
      metal <- within(metal,{P.value=format(2*pnorm(mpfr(-abs(z),100),lower.tail=TRUE,log.p=FALSE))})
      write.table(metal,file=paste0(prot,".p"),sep=\"\\t\",row.names=FALSE,quote=FALSE)
    END
    awk "NR>1" {}.p
  )'
}

function annotate()
{
  echo "--> Variant annotation"
  export annovar_home=/services/tools/annovar/2018apr16
  export humandb=$annovar_home/humandb
  export example=$annovar_home/example
  $annovar_home/annotate_variation.pl --geneanno -otherinfo -buildver hg19 $example/ex1.avinput $humandb/ --outfile ex1
}

function nodup()
# notes from https://www.biostars.org/p/264584/
{
  echo "--> remove duplicates"
  export LC_ALL=C
  (
    grep '^#' input.vcf
    grep -v "^#" input.vcf | \
    sort -t $'\t' -k1,1 -k2,2n -k4,4 | \
    awk -F '\t' 'BEGIN {prev="";} {key=sprintf("%s\t%s\t%s",$1,$2,$4);if(key==prev) next;print;prev=key;}'
  )
}

function lambda()
{
  echo "--> GC lambda"
  (
  ls METAL/*-1.tbl.gz | \
  sed 's|METAL/||g;s/-1.tbl.gz//g' | \
  parallel -j4 --env rt -C' ' '
    gunzip -c METAL/{}-1.tbl.gz | \
    cut -f3,12 | \
    gzip -f > work/INF1.{}.p.gz
    export protein={}
    R --no-save -q <<\ \ END
      library(gap)
      protein <- Sys.getenv("protein")
      gz <- gzfile(paste0("work/INF1.",protein,".p.gz"))
      p <- read.delim(gz,as.is=TRUE)
      cat(protein,"GC.lambda=",gc.lambda(with(p,P.value)),"\n")
  END'
  ) > work/INF1.lambda.log
  grep GC.lambda work/INF1.lambda.log | \
  grep -v gc.lambda | \
  sed 's/GC.lambda=//g' > work/INF1.lambda.dat
}

$1

# 2-5-2020 JHZ

export TMPDIR=/rds/user/jhz22/hpc-work/work
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export tag=_nold

for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do

echo $p
export p=${p}
(
  mergeBed -i ${INF}/sentinels/${p}_nold.p -d 1000000 -c 13 -o min | \
  awk -v OFS="\t" -v prot=${p} '
  {
    if(NR==1) print "Chrom", "Start", "End", "logP", "prot"
    print $0, prot
  }'
) > ${INF}/work/${p}.merged
(
  cut -f1-4,13 ${INF}/sentinels/${p}_nold.p | \
  bedtools intersect -a ${INF}/work/${p}.merged -b - -wa -wb | \
  awk '$4==$10' | \
  cut -f1-5,9,10 | \
  awk -v OFS="\t" '
  !(/CCL25/&&/chr19:49206145_C_G/){
    if(NR==1) print "Chrom", "Start", "End", "log10p", "prot", "MarkerName", "log10p_check", "CHR", "POS", "SNP", "BP"
    CHR=substr($1,4)
    split($6,noalleles,"_")
    split(noalleles[1],chrpos,":")
    POS=chrpos[2]
    SNP=$6
    BP=chrpos[2]
    print $0,CHR,POS,SNP,BP
  }'
) > ${INF}/work/${p}.sentinels

done

(
  cat ${INF}/work/*sentinels | head -1
  for p in $(ls sentinels/*${tag}.p | sed 's|sentinels/||g;s|'"$tag"'.p||g'); do awk 'NR>1' ${INF}/work/${p}.sentinels; done
) > INF1.merge

R --no-save -q <<END
  require(gap)
  tag <- Sys.getenv("tag")
  rt <- paste0("INF1")
  clumped <- read.table(paste0(rt,".merge"),as.is=TRUE,header=TRUE)
  hits <- merge(clumped[c("CHR","POS","MarkerName","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
  cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
  cis.vs.trans <- with(cistrans,data)
  write.table(cis.vs.trans,file=paste0(rt,".merge.cis.vs.trans"),row.names=FALSE,quote=TRUE)
  cis <- subset(cis.vs.trans,cis.trans=="cis")["SNP"]
  write.table(cis,file=paste0(rt,".merge.cis"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  sink(paste0(rt,".merge.out"))
  with(cistrans,table)
  sink()
  with(cistrans,total)
  pdf(paste0(rt,".merge.circlize.pdf"))
  circos.cis.vs.trans.plot(hits=paste0(rt,".merge"),inf1,"uniprot")
  dev.off()
END

R --no-save -q <<END
  merge.out <- read.table("INF1.merge.out",as.is=TRUE,header=TRUE)
  n <- nrow(merge.out)
  merged <- merge.out[-n,]
  cis.only <- with(merged,cis>0 & trans==0)
  trans.only <- with(merged,cis==0 & trans>0)
  both.cis.and.trans <- with(merged, cis>0 & trans>0)
  nrow(merged[cis.only,])
  nrow(merged[trans.only,])
  nrow(merged[both.cis.and.trans,])
END

pdftopng -r 300 INF1.merge.circlize.pdf INF1.merge.circlize
mv INF1.merge.circlize-000001.png INF1.merge.circlize.png

R --no-save -q <<END
  library(gap)
  d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  pdf("INF1.merge.pdf")
  mhtplot2d(d)
  dev.off()
END

pdftopng -r 300 INF1.merge.pdf INF1.merge
mv INF1.merge-000001.png INF1.merge.png

# rsid
awk 'NR>1' work/INF1.merge | cut -f6 | sort -k1,1 | uniq | \
join - work/INTERVAL.rsid > work/INF1.merge.rsid
cut -d' ' -f2 work/INF1.merge.rsid > work/INF1.merge.snp
grep -f work/INF1.merge.cis -v work/INF1.merge.rsid > work/INF1.merge.trans.rsid
cut -d' ' -f2 work/INF1.merge.trans.rsid > work/INF1.merge.trans.snp

# snpid --> rsid
for f in INF1.merge INF1.merge.cis.vs.trans
do
  cp ${f} ${f}-rsid
  (
  cat INF1.merge.rsid | \
  parallel --dry-run -C' ' "
    export s={1};
    export r={2};
    sed -i 's/'\"\${s}\"'/'\"\${r}\"'/g' ${f}-rsid
  "
  ) | bash
done

(
  bedtools intersect -a work/INF1.merge -b tryggve/high-LD-regions-hg19.bed | \
  sortBed | \
  mergeBed -i - -d 1000000 | \
  sed 's/chr//g' | \
  awk '{print 0 $1 ":" $2 "-" $3}'
) > ukb/ukb.range

sed '1d' work/INF1.merge | cut -f5 | sort -k1,1 | uniq | join -t$'\t' - work/inf1.tmp | sort -k2,2 > work/INF1.merge.prot
(
  echo -e "uniprot\tprot\ttarget\tgene"
  sed '1d' doc/olink.inf.panel.annot.tsv | \
  sed 's/\"//g' | \
  cut -f2,3,7 | \
  sort -t$'\t' -k2,2 | \
  join -j2 -t$'\t' work/INF1.merge.prot -
) > work/INF1.merge.id
sed '1d' work/INF1.merge.id | \
awk '!/NA/' | \
cut -f4 > work/INF1.merge.gene

export UKB=/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/imputed/uk10k_hrc/HRC_UK10K
(
  grep '#' -h -v $UKB/ukb_impv3_chr*_snpstats.txt | \
  head -1 | \
  awk '{$1=$1 "\t" "id"};1'
  for ((i=1;i<23;i++))
  do
    grep '#' -v  $UKB/ukb_impv3_chr${i}_snpstats.txt | \
    grep -w -f work/INF1.merge.snp 
  done
) | \
sed 's/ /\t/g' | \
awk -v OFS="\t" '{
  chr=$3+0;pos=$4;a1=$5;a2=$6
  if (a1>a2) snpid="chr" chr ":" pos "_" a2 "_" a1;
  else snpid="chr" chr ":" pos "_" a1 "_" a2
  if (NR>1) $1=snpid "\t" $1;
};1' > work/INF1.merge.ukbsnp

R --no-save -q <<END
  library(gap)
# library(Rmpfr)
  gz <- gzfile("METAL/IL.17C-1.tbl.gz")
  IL.17C <- within(read.delim(gz,as.is=TRUE), {
   Z <- Effect/StdErr;
#  P <- as.numeric(2*pnorm(mpfr(abs(Z),100),lower.tail=FALSE))
   P <- 2*pnorm(abs(Z),lower.tail=FALSE)
  })
  subset(IL.17C, P==0)
  png("IL.17C.png", res=300, units="in", width=9, height=6)
  par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
  mhtplot.trunc(IL.17C, chr="Chromosome", bp="Position", p="P", snp="MarkerName", z = "Z",
                suggestiveline=FALSE, genomewideline=-log10(5e-10), logp = TRUE,
                cex.mtext=0.6, cex.text=0.7,
                mtext.line=4, y.brk1=300, y.brk2=500, cex.axis=0.6, cex.y=0.6, cex=0.5,
                y.ax.space=20,
                col = c("blue4", "skyblue")
  )
  dev.off()
END

R --no-save <<END
  library(gap)
  d <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  r <- mhtplot2d(d)
  r <- within(r,{z=-z})
  head(r)
  write.table(r,"INF1.merge.plotly",quote=FALSE,row.names=FALSE,sep=",")
  r <- within(r,{x=x/1e9;y=y/1e9;z=z/1e2})
  write.csv(subset(r,col=="red"),"red.dat",quote=FALSE,row.names=FALSE)
  write.csv(subset(r,col=="blue"),"blue.dat",quote=FALSE,row.names=FALSE)
END
paste -d',' <(cut -d',' -f1-3 blue.dat) <(cut -d',' -f1-3 red.dat) | \
awk -v FS="," -v OFS="," '{if(NF==4) print $1,$2,$3,",,"; else print}' |
awk '{if(NR==1) print "x1,y1,z1,x2,y2,z2"; else print}' > INF1.merge.d3

awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2>0 && $3==0' | wc -l
awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2==0 && $3>0' | wc -l
awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2>0 && $3>0' | wc -l

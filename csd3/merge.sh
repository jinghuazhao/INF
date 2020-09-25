#!/usr/bin/bash

export TMPDIR=/rds/user/jhz22/hpc-work/work
export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export tag=_nold

module load gcc/6

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
  clumped <- read.delim(paste0(rt,".merge"),as.is=TRUE)
  hits <- merge(clumped[c("CHR","POS","MarkerName","prot","log10p")],inf1[c("prot","uniprot")],by="prot")
  names(hits) <- c("prot","Chr","bp","SNP","log10p","uniprot")
  cistrans <- cis.vs.trans.classification(hits,inf1,"uniprot")
  cis.vs.trans <- with(cistrans,data)
  write.csv(cis.vs.trans,file=paste0(rt,".merge.cis.vs.trans"),quote=FALSE,row.names=FALSE)
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

awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2>0 && $3==0' | wc -l
awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2==0 && $3>0' | wc -l
awk 'NR==2,NR==71' work/INF1.merge.out | awk '$2>0 && $3>0' | wc -l

pdftopng -r 300 INF1.merge.circlize.pdf INF1.merge.circlize
mv INF1.merge.circlize-000001.png INF1.merge.circlize.png

R --no-save -q <<END
  library(gap)
  d <- read.csv("INF1.merge.cis.vs.trans",as.is=TRUE)
  png("INF1.merge.png",height=20,width=20,units="cm",res=300)
  mhtplot2d(d)
  dev.off()
END

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
  d <- read.csv("INF1.merge.cis.vs.trans",as.is=TRUE)
  r <- mhtplot2d(d)
  r <- within(r,{z=-z})
  head(r)
  write.csv(r,"INF1.merge.plotly",quote=FALSE,row.names=FALSE)
  r <- within(r,{x=x/1e9;y=y/1e9;z=z/1e2})
  write.csv(subset(r,col=="red"),"red.dat",quote=FALSE,row.names=FALSE)
  write.csv(subset(r,col=="blue"),"blue.dat",quote=FALSE,row.names=FALSE)
END
paste -d',' <(cut -d',' -f1-3 blue.dat) <(cut -d',' -f1-3 red.dat) | \
awk -v FS="," -v OFS="," '{if(NF==4) print $1,$2,$3,",,"; else print}' |
awk '{if(NR==1) print "x1,y1,z1,x2,y2,z2"; else print}' > INF1.merge.d3

function INTERVAL()
# SCALLOP/INF -- INTERVAL overlap
{
  export OLINK=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/jp549/olink-merged-output
  ls $OLINK/*gz | xargs -l basename -s _chr_merged.gz | grep -v -e cvd -e P23560 | sed 's/INTERVAL_inf1_//;s/___/ /'> INTERVAL.list
  (
    gunzip -c ${OLINK}/INTERVAL_cvd3_SELP___P16109_chr_merged.gz | \
    awk 'NR==1{print "UniProt","prot","chr","pos",$2,$22,$24,$25}'
    join <(awk 'NR>1{print $5,$8 ":" $9}' work/INF1.merge | sort -k1,1) <(sort -k1,1 INTERVAL.list) | \
    awk '{
       gsub(/chr/,"",$2);
       split($2,a,":");
       chr=a[1];
       pos=a[2];
       print $0,chr,pos
    }' | \
    parallel --env OLINK -C' ' '
      zgrep -H -w {5} ${OLINK}/INTERVAL_inf1_{1}___{3}_chr_merged.gz | \
      awk -vchr={4} "(\$3==chr)" | \
      awk -vuniprot={3} -v prot={1} -v chr={4} -v pos={5} "{print uniprot, prot, chr, pos, \$2,\$22,\$24,\$25}"
    '
  ) > INTERVAL.overlap

  join -23 <(sort -k1,1 INTERVAL.list) <(awk 'NR>1{print $2,$3,$5,$6,$8 ":" $9}' work/INF1.merge | sort -k3,3) | \
  awk '{
         if($3==$4) {$3=$3-1e6;$4=$4+1e6}
         if($3<0) $3=0
         gsub(/chr/,"",$6);
         split($6,a,":");
         chr=a[1];
         print $0,chr
  }' | \
  parallel -j5 --env OLINK -C' ' '
  (
     gunzip -c ${OLINK}/INTERVAL_cvd3_SELP___P16109_chr_merged.gz | \
     awk "NR==1{print \"UniProt\",\"prot\",\"rsid\",\"chr\",\"pos\",\$22,\$24,\$25}"
     zcat ${OLINK}/INTERVAL_inf1_{1}___{2}_chr_merged.gz | \
     awk -v chr={7} -v start={3} -v end={4} -v NA="NA" "chr==\$3+0 && \$4>=start && \$4<=end && index(\$0,NA)==0" | \
     awk -v prot={1} -v uniprot={2} "{print uniprot, prot, \$2, \$3+0, \$4,\$22,\$24,\$25}"
  ) | \
   gzip -f > INTERVAL.merge.{2}-{1}-{5}.gz
  '
}

INTERVAL

# region according to INF1
join -11 -25 <(sort -k1,1 work/inf1.tmp | grep -v P23560) \
             <(sed '1d' work/INF1.merge | sort -k5,5 | \
               awk '{if($3-$2==1) {$2=$2-1e6;$3=$3+1e6};
                     if($2<0) $2=0;
                     print}') | \
cut -d' ' -f1-2,7,9,10 | \
parallel -j5 -C' ' '
  gunzip -c METAL/{1}-1.tbl.gz | \
  awk -vchr={4} -vpos={5} "NR==1||(\$1==chr && \$2>=pos-1e6 && \$2<pos+1e6)" | \
  cut -f1-6,10-12 | \
  gzip -f > work/INF1.merge.{2}-{1}-{3}.gz
'

# regions according to SomaLogic
R --no-save -q <<END
  library(pQTLtools)
  sentinels <- st4[,5:12]
  snpid <- gap::chr_pos_a1_a2(st4[,7],st4[,8],st4[,11],st4[,12])
  write.table(cbind(snpid,sentinels),file="SomaLogic.sentinels",quote=FALSE,row.names=FALSE)
END

# --- protein overlap
# Somalogic proteins with sentinels (1469) - NOTE P29460,Q9NPF7 in SomaLogic
sed '1d' SomaLogic.sentinels | sort -k2,2 | cut -d' ' -f2 | uniq | wc -l
cut -f2 work/inf1.tmp | grep -v P23560 > work/INF1.uniprot

# number of proteins with sentinels in both Olink and SomaLogic (28)
cut -f2 work/INF1.merge.prot | grep -f - SomaLogic.sentinels | cut -d' ' -f2 | sort | uniq | wc -l
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - work/INF1.merge.prot | wc -l

# --- signal overlap
# all SomaLogic signals in Olink (51)
join -j2 <(sort -k2,2 work/inf1.tmp | grep -v P23560) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) | wc -l
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f work/INF1.uniprot - | wc -l

# all SomaLogic signals from overlapping proteins (45)
sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2 | grep -f work/INF1.merge.uniprot - | wc -l
join -j2 <(sort -k2,2 work/INF1.merge.prot) <(sed '1d;s/,Q9NPF7//' SomaLogic.sentinels | sort -k2,2) > SomaLogic.INF1.all
(
cut -d' ' -f1-3,5,7,8 SomaLogic.INF1.all | \
parallel -C' ' '
  zgrep -w {3} METAL/{2}-1.tbl.gz
  gunzip -c METAL/{2}-1.tbl.gz | \
  awk -vchr={4} -vstart={5} -vend={6} "NR==1||(\$1==chr&&\$2>=start&&\$2<=end&&\$12<-9.30103)" | \
  cut -f1-5,10-12 | \
  gzip -f > INF1.SomaLogic.{1}-{2}-{3}.gz
  export lines=$(gunzip -c INF1.SomaLogic.{1}-{2}-{3}.gz | wc -l | cut -d" " -f1)
  if [ ${lines} -eq 1 ]; then rm INF1.SomaLogic.{1}-{2}-{3}.gz; fi
'
) > INF1.SomaLogic.all

# which are also genomewide significant (32)
awk '$12<-9.30103' INF1.SomaLogic.all | wc -l

# identical signals (10)
cat SomaLogic.INF1.all | \
parallel -C' ' 'awk -v prot={2} -v MarkerName={3} "\$5==prot && \$6==MarkerName" work/INF1.merge'
join <(awk '{print $2"-"$3,$0}' SomaLogic.INF1.all | sort -k1,1) <(awk '{print $5"-"$6,$0}' work/INF1.merge | sort -k1,1)

# Olink overlapping proteins
cut -d' ' -f2 SomaLogic.sentinels | sed 's/P29460,Q9NPF7/P29460/' | grep -f - work/INF1.merge.prot | \
cut -f1 | grep -f - work/INF1.merge | \
cut -f5 | sort | uniq | wc -l

R --no-save -q <<END
  library(VennDiagram)
  INF1_prot <- read.table("work/INF1.merge.prot",col.names=c("prot","uniprot"))
  library(pQTLtools)
  SomaLogic_prot <- unique(replace(st4$UniProt,st4$UniProt=="P29460,Q9NPF7","P29460"))
  plist <- list(INF1_prot$uniprot,SomaLogic_prot)
  ov <- VennDiagram::calculate.overlap(plist)
  ov$a3
  venn.plot <- draw.pairwise.venn(1469, 70, 28,
               category = c("SomaLogic", "Olink"),
               fill = c("blue", "red"),
               lty = "blank",
               cex = 2,
               cat.cex = 2,
               cat.pos = c(200, 50),
               cat.dist = 0.09,
               cat.just = list(c(-1, -1), c(1, 1)),
               ext.pos = 30,
               ext.dist = -0.05,
               ext.length = 0.85,
               ext.line.lwd = 2,
               ext.line.lty = "dashed"
             );
  grid.draw(venn.plot);
  png('SomaLogic-Olink-proteins.png', height=20, width=20, units="cm", res=300)
  grid.draw(venn.plot);
  dev.off();
  grid.newpage()
  venn.plot <- draw.pairwise.venn(45, 83, 10,
               category = c("SomaLogic", "Olink"),
               fill = c("blue", "red"),
               lty = "blank",
               cex = 2,
               cat.cex = 2,
               cat.pos = c(30, 10),
               cat.dist = 0.09,
               cat.just = list(c(-1, -1), c(1, 1)),
               ext.pos = 30,
               ext.dist = -0.05,
               ext.length = 0.85,
               ext.line.lwd = 2,
               ext.line.lty = "dashed"
             );
  grid.draw(venn.plot);
  png('SomaLogic-Olink-sentinels.png', height=20, width=20, units="cm", res=300)
  grid.draw(venn.plot);
  dev.off();
  grid.newpage()
  st6[c(39,53),"UniProt"] <- "P29460"
  st6ov <- subset(st6,UniProt%in%ov$a3)
# significant on INTERVAL (27)
  dim(subset(st6ov,as.numeric(p.1)<1e-11))
  z <- with(st6ov,{
    chr <- st6ov[["Chr"]]
    pos <- st6ov[["Pos"]]
    a1 <- st6ov[["Effect.Allele.(EA)"]]
    a2 <- st6ov[["Other.Allele.(OA)"]]
    cbind(UniProt,snpid=gap::chr_pos_a1_a2(chr,pos,a1,a2))
  })
  write.table(merge(inf1,z,by.x="uniprot",by.y="UniProt")[c("snpid","prot","uniprot")],
              file="SomaLogic.id3",col.names=FALSE,row.names=FALSE,quote=FALSE)
END

# --- SomaLogic --> Olink lookup --- NOTE these were not necessary Olink sentinels
# Similar to ST6 of the SomaLogic paper and pQTLtools/inst/scripts/STs.R
# from the replicates how many were from INF1 (41)
(
cat SomaLogic.id3 | \
parallel -C' ' '
  zgrep -w {1} METAL/{2}-1.tbl.gz
'
) > SomaLogic.INF1
# how many among INF1 variants were genomewide significant (28)
wc -l SomaLogic.INF1
awk '$12<-9.30103' SomaLogic.INF1 | wc -l

(
cat SomaLogic.id3 | \
parallel -C' ' '
  zgrep -w {1} METAL/{2}-1.tbl.gz | \
  awk -v prot={2} -v uniprot={3} "{print prot,uniprot,\$0}"
'
) | \
sort -k5,5 | \
join -a1 -15 -e "NA" - work/INTERVAL.rsid > SomaLogic.INF1-rsid
awk '$14<-9.30103 {print $2, $21}' SomaLogic.INF1-rsid

rm SomaLogic.id3 SomaLogic.INF1 INF1.SomaLogic*gz SomaLogic.INF1.all  SomaLogic.INF1-rsid  SomaLogic.sentinels INF1.SomaLogic.all

# REACTOME
cut -d, -f10,14 work/INF1.merge.cis.vs.trans | \
sed '1d' | \
grep cis | \
cut -d, -f1 | \
sort | \
uniq | \
xsel -i

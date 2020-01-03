# 3-1-2020 JHZ

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

cd work
cut -f8,9,10 INF1.merge | \
awk -vOFS="\t" 'NR>1{split($3,a,"_");print $1,$2,$2,a[2],a[3]}' | \
sort -k1,1n -k2,2n | \
uniq > INF1.merge.avinput
R --no-save -q <<END
  cvt <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  ord <- with(cvt,order(Chr,bp))
  trans <- subset(cvt[ord,],cis.trans=="trans")
  s <- with(trans,unique(gap::inv_chr_pos_a1_a2(SNP,prefix="")))
  vars <- c("chr","pos","pos","a1","a2")
  write.table(s[vars],file="INF1.merge.trans.avinput",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
  vepinput <- "INF1.merge.trans.vepinput"
  cat("##fileformat=VCFv4.0\n", file=vepinput)
  cat("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO\n",file=vepinput,append=TRUE,sep="\t")
  s <- within(s,{snp <- paste0(chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
  vars <- c("chr","pos","snp","a1","a2","qual","filter","info")
  write.table(s[vars],file=vepinput,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
END
vep -i INF1.merge.trans.vepinput -o INF1.merge.trans.vepoutput --force_overwrite --offline

(
  cut -f1,2 ${INF}/work/INF1.merge.trans.vepinput | \
  awk -v OFS="\t" -v flanking=500000 'NR>2 {
    if ($2-flanking<0) print $1, 0, $2+flanking;
    else print $1, $2-flanking, $2+flanking
  }'
) > a1

(
  sort -k1,1n -k2,2n ${INF}/csd3/glist-hg19 | \
  grep -v X | \
  grep -v Y | \
  awk '{$1="chr" $1;print}' | \
  sed 's/ /\t/g'
) > a2

bedtools intersect -a a1 -b a2 -wa -wb -loj | \
cut  -f1-3,7 > INF1.merge.glist-hg19

R --no-save -q <<END
  vo <-  read.delim("INF1.merge.trans.vepoutput",skip=40)
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
  attr <- listAttributes(ensembl)
  filter <- listFilters(ensembl)
  s1 <- c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'description', 'hgnc_symbol', 'entrezgene_id')
  s2 <- c('ensembl_peptide_id', 'peptide', 'protein_id')
  s3 <- c('clinical_significance')
  gene <- getBM(attributes = s1, mart = ensembl)
  vepbiomart <- merge(vo,gene,by.x="Gene",by.y="ensembl_gene_id")
  write.table(vepbiomart,file="INF1.merge.trans.vepbiomart",row.name=FALSE,quote=FALSE,sep="\t")
  library(openxlsx)
  xlsx <- "INF1.merge.trans.vepbiomart.xlsx"
  unlink(xlsx, recursive = FALSE, force = TRUE)
  wb <- createWorkbook(xlsx)
  snpid_rsid <- read.table("INF1.merge.rsid",col.names=c("snpid","rsid"))
  d <- merge(snpid_rsid,vepbiomart,by.x="snpid",by.y="X.Uploaded_variation",all.y=TRUE)
  addWorksheet(wb, "vepbiomart")
  writeDataTable(wb, "vepbiomart",d)
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
END
export humandb=$HPC_WORK/annovar/humandb
# on tryggve
module load perl/5.24.0 annovar/2019oct24
export annovar_home=/services/tools/annovar/2019oct24
export humandb=$annovar_home/humandb
for s in INF1.merge INF1.merge.trans
do
  annotate_variation.pl --genomebinsize 500k -buildver hg19 ${s}.avinput $humandb/ --outfile ${s}
  convert2annovar.pl -format annovar2vcf ${s}.avinput > ${s}.vcf
  vep -i ${s}.vcf -o ${s}.vcfoutput --offline
done
cd -

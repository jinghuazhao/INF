# 12-3-2020 JHZ

export bfile=${INF}/INTERVAL/cardio/INTERVAL

cd work
cut -f1-3,5,6,9 INF1.merge | \
awk -v d=10000000 'NR>1 {gsub(/chr/,"",$1);if($6>d) from=$6-d; else from=0; to=$6+d;print $1,from,to,$4,$5}' | \
parallel --env bfile=$bfile -j5 -C' ' '
plink --bfile ${bfile} --ld-snp {5} --chr {1} --from-bp {2} --to-bp {3} --r2 --ld-window-r2 0.6 --out {4}-{5}
'
# all variants
# awk -v chr={1} -v from={2} -v to={3} "chr==\$1 && \$4 >= from && \$4 <= to" ${bfile}.bim

(
  echo prot SNPID CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
  cut -f1-3,5,6,9 --output-delimiter=' ' INF1.merge | \
  sed "1d" | \
  parallel -C' ' 'awk -v prot={4} -v snpid={5} "NR>1 && !(\$3==\$6) {print prot,snpid,\$1,\$2,\$3,\$4,\$5,\$6,\$7}" {4}-{5}.ld'
) > INF1.proxy.ld

R --no-save -q <<END
  ld <- read.table("INF1.proxy.ld",as.is=TRUE,header=TRUE)
  ord <- with(ld,order(CHR_B,BP_B))
  ct <- ld[ord,]
  all <- with(ct,unique(gap::inv_chr_pos_a1_a2(SNP_B,prefix="")))
  all <- within(all,{snp <- paste0(chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
  writetable <- function(d,f,...) write.table(d,file=f,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",...)
  for (f in c("INF1.proxy"))
  {
    avinput <- paste0(f,".avinput")
    vars <- c("chr","pos","pos","a1","a2")
    writetable(all[vars],avinput)
    vepinput <- paste0(f,".vepinput")
    cat("##fileformat=VCFv4.0\n", file=vepinput)
    cat("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO\n",file=vepinput,append=TRUE,sep="\t")
    vars <- c("chr","pos","snp","a1","a2","qual","filter","info")
    writetable(all[vars],vepinput,append=TRUE)
  }
END

export ANNOVAR=${HPC_WORK}/annovar
annotate_variation.pl -buildver hg19 INF1.proxy.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile INF1.proxy
vep -i INF1.proxy.vepinput -o INF1.proxy.vepoutput --pick --check_existing --distance 500000 --force_overwrite --offline --everything --assembly GRCh37

cd -

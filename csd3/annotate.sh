# 23-1-2020 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export ANNOVAR=${HPC_WORK}/annovar
export POLYPHEN=$HPC_WORK/polyphen-2.2.2
export VEP=${HPC_WORK}/ensembl-vep
export TMPDIR=/rds/user/jhz22/hpc-work/work

cd work
R --no-save -q <<END
  cvt <- read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  ord <- with(cvt,order(Chr,bp))
  ct <- cvt[ord,]
  all <- with(ct,unique(gap::inv_chr_pos_a1_a2(SNP,prefix="")))
  all <- within(all,{snp <- paste0(chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
  trans <- with(subset(ct,cis.trans=="trans"),unique(gap::inv_chr_pos_a1_a2(SNP,prefix=""))) 
  trans <- within(trans,{snp <- paste0(chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
  writetable <- function(d,f,...) write.table(d,file=f,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",...)
  for (f in c("INF1.merge","INF1.merge.trans"))
  {
    avinput <- paste0(f,".avinput")
    vars <- c("chr","pos","pos","a1","a2")
    if(f=="INF1.merge") writetable(all[vars],avinput) else writetable(trans[vars],avinput)
    vepinput <- paste0(f,".vepinput")
    cat("##fileformat=VCFv4.0\n", file=vepinput)
    cat("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO\n",file=vepinput,append=TRUE,sep="\t")
    vars <- c("chr","pos","snp","a1","a2","qual","filter","info")
    if(f=="INF1.merge") writetable(all[vars],vepinput,append=TRUE) else writetable(trans[vars],vepinput,append=TRUE)
  }
END
(
  head -1 INF1.merge
  cut -f6  INF1.merge | sed '1d' | sort -k1,1 | uniq | grep -f INF1.merge.cis -v > INF1.merge.notcis
  grep -f INF1.merge.notcis INF1.merge
) > INF1.merge.trans
for s in INF1.merge INF1.merge.trans
do
   export s=${s}
   # ANNOVAR
   annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
   table_annovar.pl ${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
        -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
        -operation g,g,g,g,r,f,f,f,r \
        -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
   # Polyphen-2
   grep -v -w -f INF1.merge.cis ${s} | \
   cut -f6 | \
   sed '1d;s/_/ /;s/_/\//' | \
   sort -k1,1 | \
   uniq > ${s}.pph.list
   mapsnps.pl -g hg19 -m -U -y ${s}.pph.input ${s}.pph.list 1>${s}.pph.features 2>${s}.log
   run_pph.pl ${s}.pph.input 1>${s}.pph.output 2>${s}.pph.log
   run_weka.pl ${s}.pph.output >${s}.pph.humdiv.output
   run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model ${s}.pph.output >${s}.pph.humvar.output
   # VEP
   vep -i ${s}.vepinput -o ${s}.vepoutput --pick --check_existing --distance 500000 --force_overwrite --offline --everything --assembly GRCh37
   vep -i ${s}.vepinput --species homo_sapiens -o ${s}.clinvar \
       --cache --distance 500000 --offline --force_overwrite \
       --assembly GRCh37 --pick --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,DBVARID,MC,RS
   vep -i ${s}.vepinput -o ${s}.dbNSFP --cache --distance 500000 --force --offline --pick \
       --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,clinvar_id,clinvar_clnsig,clinvar_review,clinvar_trait,1000Gp3_EUR_AF,CADD_phred,Eigen-PC-phred_coding,ExAC_NFE_AF,LRT_pred,FATHMM_pred,GERP++_RS,GTEx_V7_tissue,MutPred_protID,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,SIFT4G_pred,fathmm-MKL_coding_pred,rs_dbSNP151,fathmm-MKL_coding_pred,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF
   vep --af_1kg --af_esp --af_gnomad --appris --biotype --buffer_size 500 \
       --ccds --check_existing --distance 500000 --domains --hgvs --mane --pick \
       --polyphen b --protein --pubmed --regulatory --sift b --species homo_sapiens \
       --symbol --transcript_version --tsl --uniprot --cache --input_file ${s}.vepinput \
       --output_file ${s}.vcf --port 3337 --vcf
   R --no-save <<\ \ \ END
     library(ensemblVEP)
     s <- Sys.getenv("s")
     f <- paste0(s,".vcf")
     vcf <- readVcf(f, "hg19")
     csq <- parseCSQToGRanges(f, VCFRowID=rownames(vcf))
     write.table(mcols(csq),file=paste0(s,".txt"), quote=FALSE, sep="\t")
   END
   vep --af_1kg --af_esp --af_gnomad --appris --biotype --buffer_size 500 \
       --ccds --check_existing --distance 500000 --domains --hgvs --mane --pick \
       --polyphen b --protein --pubmed --regulatory --sift b --species homo_sapiens \
       --symbol --transcript_version --tsl --uniprot --cache --input_file ${s}.vepinput \
       --output_file ${s}.vepweb --port 3337
done

export skips=$(grep '##' INF1.merge.trans.vepoutput | wc -l)
R --no-save -q <<END
  cvt <- subset(read.table("INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE),cis.trans=="trans")
  vf <- read.delim("INF1.merge.trans.variant_function",header=FALSE)
  names(vf) <- c("function","ensGene","chr","pos","pos2","a1","a2")
  evf <- read.delim("INF1.merge.trans.exonic_variant_function",header=FALSE)
  names(evf) <- c("lineno","substitution","description","chr","pos","pos2","a1","a2")
  vf <- within(vf,{SNP <- paste0(chr,":",pos,"_",a1,"_",a2)})
  evf <- within(evf,{SNP <- paste0(chr,":",pos,"_",a1,"_",a2)})
  cvtvf <- merge(cvt,vf[c("function","ensGene","SNP")],by="SNP",all=TRUE)
  cvtvfevf <- merge(cvtvf,evf[c("substitution","description","SNP")],by="SNP",all=TRUE)
  vo <-  read.delim("INF1.merge.trans.vepoutput",skip=as.integer(Sys.getenv("skips")))
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
  attr <- listAttributes(ensembl)
  filter <- listFilters(ensembl)
  s1 <- c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'description', 'hgnc_symbol', 'entrezgene_id')
  s2 <- c('ensembl_transcript_id', 'transcription_start_site', 'transcript_start', 'transcript_end')
  s3 <- c('ensembl_peptide_id', 'peptide', 'protein_id')
  s4 <- c('clinical_significance')
  gt <- getBM(attributes = s1, mart = ensembl)
  vepbiomart <- merge(vo,gt,by.x="Gene",by.y="ensembl_gene_id")
  write.table(vepbiomart,file="INF1.merge.trans.vepbiomart",row.name=FALSE,quote=FALSE,sep="\t")
  library(openxlsx)
  xlsx <- "INF1.merge.trans.annovarvepbiomart.xlsx"
  unlink(xlsx, recursive = FALSE, force = TRUE)
  wb <- createWorkbook(xlsx)
  snpid_rsid <- read.table("INF1.merge.rsid",col.names=c("snpid","rsid"))
  d <- merge(snpid_rsid,vepbiomart,by.x="snpid",by.y="X.Uploaded_variation",all.y=TRUE)
  addWorksheet(wb, "annovar")
  writeDataTable(wb, "annovar", cvtvfevf)
  addWorksheet(wb, "vepbiomart")
  writeDataTable(wb, "vepbiomart",d)
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
END

# sentinel positions +/- 500k
(
  cut -f1,2 ${INF}/work/INF1.merge.trans.vepinput | \
  awk -v OFS="\t" -v flanking=500000 'NR>2 {
    if ($2-flanking<0) print "chr" $1, 0, $2+flanking;
    else print "chr" $1, $2-flanking, $2+flanking
  }'
) > a1
(
  sort -k1,1n -k2,2n ${INF}/csd3/glist-hg19 | \
  grep -v -e X -e Y -w | \
  awk '{$1="chr" $1;print}' | \
  sed 's/ /\t/g'
) > a2
bedtools intersect -a a1 -b a2 -wa -wb -loj | \
cut  -f1-3,7 > INF1.merge.trans.glist-hg19
rm a1 a2

# ProGeM
# Bottom-up
# https://www.gtexportal.org/home/datasets
# Top-down (GenomicRanges, KEGGREST, STRINGdb)
# ftp://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
# http://www.omim.org/
# https://www.humanmine.org/
# http://www.genome.jp/kegg/
# http://string-db.org/
# BiocManager::install("garfield")

cd -

# 20-1-2020 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export ANNOVAR=${HPC_WORK}/annovar
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

for s in INF1.merge INF1.merge.trans
do
   annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
   vep -i ${s}.vepinput -o ${s}.vepoutput --pick --force_overwrite --offline --everything --assembly GRCh37
   vep -i ${s}.vepinput --species homo_sapiens -o ${s}.clinvar \
       --cache --offline --force_overwrite \
       --assembly GRCh37 --pick --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
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
    if ($2-flanking<0) print $1, 0, $2+flanking;
    else print $1, $2-flanking, $2+flanking
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

cd -

# Bottom-up
# https://www.gtexportal.org/home/datasets
# Top-down (GenomicRanges, KEGGREST, STRINGdb)
# ftp://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
# http://www.omim.org/
# https://www.humanmine.org/
# http://www.genome.jp/kegg/
# http://string-db.org/
# BiocManager::install("garfield")

# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
R --no-save <<END
  library(rentrez)
  library(XML)
# set_entrez_key("")
  Sys.getenv("ENTREZ_KEY")
  entrez_dbs()
  entrez_db_links("pubmed")
  pubmed_fields <- entrez_db_searchable("pubmed")
  r <- entrez_search(db="pubmed",term="(pQTLs OR (protein quantitative locus)) AND 2008:2020[DP] AND homo_sapiens[ORGN]",use_history=TRUE)
  class(r)
  names(r)
  with(r,web_history)
  unlink(paste("pubmed",c("fetch","summary"),sep="."))
  for(i in seq(1,with(r,count),50))
  {
    cat(i+49, "records downloaded\r")
    f <- entrez_fetch(db="pubmed", web_history=with(r,web_history), rettype="text", retmax=50, retstart=i)
    write.table(f, col.names=FALSE, row.names=FALSE, file="pubmed.fetch", append=TRUE)
    s <- entrez_summary(db="pubmed", web_history=with(r,web_history), rettype="text", retmax=50, retstart=i)
    fields <- c("uid", "pubdate", "sortfirstauthor", "title", "source", "volume", "pages")
    e <- extract_from_esummary(s, fields)
    write.table(t(e), col.names=FALSE, row.names=FALSE, file="pubmed.summary", append=TRUE, sep="\t")
  }
  link.example <- function()
  {
    upload <- entrez_post(db="omim", id=600807)
    asthma_variants <- entrez_link(dbfrom="omim", db="clinvar", cmd="neighbor_history", web_history=upload)
    asthma_variants
    snp_links <- entrez_link(dbfrom="clinvar", db="snp",
                             web_history=asthma_variants$web_histories$omim_clinvar,
                             cmd="neighbor_history")
    all_links <- entrez_link(dbfrom='pubmed', id=351, db='all')
  }
END

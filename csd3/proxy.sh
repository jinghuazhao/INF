#!/usr/bin/bash

export bfile=${INF}/INTERVAL/cardio/INTERVAL
export TMPDIR=${HPC_WORK}/work

# replace INF1.merge.trans with INF1.merge for all sentinels, add grep -f INF1.merge.cis | \ for cis sentinels.

cd work
awk 'NR>1{print $1,$5,$6,$9}' INF1.merge.trans | \
parallel --env bfile --env TMPDIR -j1 -C' ' '
  awk -vchr={1} -vpos={4} -vd=500000 "chr==\"chr\"\$1 && \$4>=pos-d && \$4<=pos+d {print \$2}" ${bfile}.bim > {2}-{3}.snplist
  plink --bfile ${bfile} --extract {2}-{3}.snplist --r2 inter-chr --ld-window-r2 0.8 --out {2}-{3}
'
(
  echo prot SNPID CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
  cut -f5,6 --output-delimiter=' ' INF1.merge.trans | \
  sed "1d" | \
  parallel -C' ' 'awk -v prot={1} -v snpid={2} "NR>1 {print prot,snpid,\$1,\$2,\$3,\$4,\$5,\$6,\$7}" {1}-{2}.ld'
) > INF1.proxy.ld

cut -f5,6 --output-delimiter=' ' INF1.merge.trans | \
sed "1d" | \
parallel -C' ' '
export ldfile={1}-{2}
R --no-save -q <<END
  ldfile <- Sys.getenv("ldfile")
  ld <- read.table(paste0(ldfile,".ld"),as.is=TRUE,header=TRUE)
  sentinels <- as.matrix(ld[c("CHR_A","BP_A","SNP_A")])
  proxies <- as.matrix(ld[c("CHR_B","BP_B","SNP_B")])
  all <- as.data.frame(rbind(sentinels,proxies))
  ord <- with(all,order(CHR_A,BP_A))
  all <- all[ord,]
  all <- with(all,unique(gap::inv_chr_pos_a1_a2(SNP_A,prefix="")))
  all <- within(all,{snp <- paste0("chr",chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
  writetable <- function(d,f,...) write.table(d,file=f,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",...)
  for (f in ldfile)
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
'

export ANNOVAR=${HPC_WORK}/annovar
export LOFTEE=${HPC_WORK}/loftee
export POLYPHEN=$HPC_WORK/polyphen-2.2.2
export VEP=${HPC_WORK}/ensembl-vep

cut -f5,6 --output-delimiter=' ' INF1.merge.trans | \
sed "1d" | \
parallel -C' ' '
export s={1}-{2}
# ANNOVAR
annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
table_annovar.pl ${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
     -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
     -operation g,g,g,g,r,f,f,f,r \
     -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
# Polyphen-2
cut -d" " -f8 ${s}.ld | \
sed "1d;s/_/ /;s/_/\//" | \
sort -k1,1 | \
uniq > ${s}.pph.list
mapsnps.pl -g hg19 -m -U -y ${s}.pph.input ${s}.pph.list 1>${s}.pph.features 2>${s}.log
run_pph.pl ${s}.pph.input 1>${s}.pph.output 2>${s}.pph.log
run_weka.pl ${s}.pph.output >${s}.pph.humdiv.output
run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model ${s}.pph.output >${s}.pph.humvar.output
# VEP
vep -i ${s}.vepinput -o ${s}.vepoutput --pick --check_existing --force_overwrite --offline --everything --assembly GRCh37 \
    --nearest symbol --symbol --pubmed --uniprot --protein --sift b --polyphen b --tab
vep -i ${s}.vepinput -o ${s}.clinvar --species homo_sapiens \
    --cache --offline --force_overwrite \
    --assembly GRCh37 --pick --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,DBVARID,MC,RS \
    --fields Uploaded_variation,Gene,Consequence,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar_DBVARID,ClinVar_MC,ClinVar_RS --tab
vep --af_1kg --af_esp --af_gnomad --appris --biotype --buffer_size 500 \
    --ccds --check_existing --domains --hgvs --mane --pick \
    --polyphen b --protein --pubmed --regulatory --sift b --species homo_sapiens \
    --symbol --transcript_version --tsl --uniprot --cache --input_file ${s}.vepinput \
    --output_file ${s}.vepweb --port 3337 --tab --force_overwrite
cd $LOFTEE
export dbNSFP_1=clinvar_id,clinvar_clnsig,clinvar_review,clinvar_trait,1000Gp3_EUR_AF,CADD_phred,Eigen-PC-phred_coding,ExAC_NFE_AF,LRT_pred,
export dbNSFP_2=FATHMM_pred,GERP++_RS,GTEx_V7_tissue,MutPred_protID,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,SIFT4G_pred,fathmm-MKL_coding_pred,
export dbNSFP_3=rs_dbSNP151,fathmm-MKL_coding_pred,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF
export dbNSFP_fields=${dbNSFP_1}${dbNSFP_2}${dbNSFP_3}
vep -i ${INF}/work/${s}.vepinput -o ${INF}/work/${s}.dbNSFP --cache --force --offline --pick --tab \
    --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz \
    --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,${dbNSFP_fields}

# cd ~/rds/public_databases/dbNSFP4.1a
# awk -vOFS="\t" '{if(NR>1) $3="."; print}'  ~/INF/work/INF1.proxy.trans.vepinput > ~/INF1.proxy.trans.vcf
# java -jar search_dbNSFP41a.jar -i ~/INF1.vcf -o ~/INF/work/INF1.proxy.trans.dbNSFP41 -v hg19
# cd -

# NB all files are moved into relevant directories
mv \
INF1.proxy.ld \
INF1.proxy.vepinput \
INF1.proxy.avinput \
INF1.proxy.variant_function \
INF1.proxy.exonic_variant_function \
INF1.proxy.pph.list \
INF1.proxy.hg19_multianno.csv \
INF1.proxy.pph.input \
INF1.proxy.pph.features \
INF1.proxy.log \
INF1.proxy.pph.log \
INF1.proxy.pph.output \
INF1.proxy.pph.humvar.output \
INF1.proxy.pph.humdiv.output \
INF1.proxy.vepoutput_summary.html \
INF1.proxy.vepoutput \
INF1.proxy.clinvar_summary.html \
INF1.proxy.clinvar \
INF1.proxy.vepweb_summary.html \
INF1.proxy.vepweb \
INF1.proxy.dbNSFP_summary.html \
INF1.proxy.dbNSFP annotate/trans

# ~/INF/work/annotate/trans/INF1.proxy.vepoutput
export skips=$(grep '##' ${INF}/work/annotate/trans/INF1.proxy.vepoutput | wc -l)
R --no-save -q <<END
  INF <- Sys.getenv("INF")
  cvt <- subset(read.csv(file.path(INF,"work","INF1.merge.cis.vs.trans"),as.is=TRUE),cis.trans=="trans")
  vo <- read.delim(file.path(INF,"work","annotate","trans","INF1.proxy.vepoutput"),skip=as.integer(Sys.getenv("skips")))
  cvt_vo <- merge(cvt,vo[c("X.Uploaded_variation","SWISSPROT","Gene","NEAREST")],by.x="SNP",by.y="X.Uploaded_variation")
  trans_anno <-merge(cvt_vo,biomaRt,by.x="NEAREST",by.y="hgnc_symbol",all.x=TRUE)[c("SNP","NEAREST","SWISSPROT","uniprotswissprot")]
  no_uniprot <- with(trans_anno,SWISSPROT=="-")
  trans_anno[no_uniprot,"SWISSPROT"] <- with(trans_anno,uniprotswissprot)[no_uniprot]
  pp <- unique(subset(trans_anno,select=-uniprotswissprot))
  save(pp,file=file.path(INF,"work","INF1.merge.trans.anno.rda"))
END

cd -

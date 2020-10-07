#!/usr/bin/bash

export bfile=${INF}/INTERVAL/cardio/INTERVAL

cd work
cut -f1-3,5,6,9 INF1.merge | \
awk -v d=500000 'NR>1 {gsub(/chr/,"",$1);if($6>d) from=$6-d; else from=0; to=$6+d;print $1,from,to,$4,$5}' | \
parallel --env bfile=$bfile -j5 -C' ' '
  plink --bfile ${bfile} --ld-snp {5} --chr {1} --from-bp {2} --to-bp {3} --r2 --ld-window-r2 0.6 --out {4}-{5}
# all variants
# awk -v chr={1} -v from={2} -v to={3} "chr==\$1 && \$4 >= from && \$4 <= to" ${bfile}.bim
'

(
  echo prot SNPID CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
  cut -f1-3,5,6,9 --output-delimiter=' ' INF1.merge | \
  sed "1d" | \
  parallel -C' ' 'awk -v prot={4} -v snpid={5} "NR>1 {print prot,snpid,\$1,\$2,\$3,\$4,\$5,\$6,\$7}" {4}-{5}.ld'
) > INF1.proxy.ld

R --no-save -q <<END
  ld <- read.table("INF1.proxy.ld",as.is=TRUE,header=TRUE)
  ord <- with(ld,order(CHR_B,BP_B))
  ct <- ld[ord,]
  all <- with(ct,unique(gap::inv_chr_pos_a1_a2(SNP_B,prefix="")))
  all <- within(all,{snp <- paste0("chr",chr,":",pos,"_",a1,"_",a2); qual <- "."; filter <- "."; info <- "."})
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
export LOFTEE=${HPC_WORK}/loftee
export POLYPHEN=$HPC_WORK/polyphen-2.2.2
export VEP=${HPC_WORK}/ensembl-vep

for s in INF1.proxy
do
   export s=${s}
   # ANNOVAR
   annotate_variation.pl -buildver hg19 ${s}.avinput ${ANNOVAR}/humandb/ -dbtype ensGene --outfile ${s}
   table_annovar.pl ${s}.avinput $ANNOVAR/test -buildver hg19 -out ${s} \
        -protocol ensGene,refGene,ccdsGene,wgEncodeGencodeBasicV19,cytoBand,exac03,avsnp147,dbnsfp30a,gwasCatalog \
        -operation g,g,g,g,r,f,f,f,r \
        -remove -nastring . -csvout -polish -xref $ANNOVA/example/gene_xref.txt
   # Polyphen-2
   cut -d' ' -f8 ${s}.ld | \
   sed '1d;s/_/ /;s/_/\//' | \
   sort -k1,1 | \
   uniq > ${s}.pph.list
   mapsnps.pl -g hg19 -m -U -y ${s}.pph.input ${s}.pph.list 1>${s}.pph.features 2>${s}.log
   run_pph.pl ${s}.pph.input 1>${s}.pph.output 2>${s}.pph.log
   run_weka.pl ${s}.pph.output >${s}.pph.humdiv.output
   run_weka.pl -l $POLYPHEN/models/HumVar.UniRef100.NBd.f11.model ${s}.pph.output >${s}.pph.humvar.output
   # VEP
   vep -i ${s}.vepinput -o ${s}.vepoutput --pick --check_existing --distance 500000 --force_overwrite --offline --everything --assembly GRCh37 \
       --nearest symbol --symbol --pubmed --uniprot --protein --sift b --polyphen b --tab
   vep -i ${s}.vepinput -o ${s}.clinvar --species homo_sapiens \
       --cache --distance 500000 --offline --force_overwrite \
       --assembly GRCh37 --pick --custom clinvar_GRCh37.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,DBVARID,MC,RS \
       --fields Uploaded_variation,Gene,Consequence,ClinVar_CLNSIG,ClinVar_CLNREVSTAT,ClinVar_CLNDN,ClinVar_DBVARID,ClinVar_MC,ClinVar_RS --tab
   vep --af_1kg --af_esp --af_gnomad --appris --biotype --buffer_size 500 \
       --ccds --check_existing --distance 500000 --domains --hgvs --mane --pick \
       --polyphen b --protein --pubmed --regulatory --sift b --species homo_sapiens \
       --symbol --transcript_version --tsl --uniprot --cache --input_file ${s}.vepinput \
       --output_file ${s}.vepweb --port 3337 --tab --force_overwrite
   cd $LOFTEE
   export dbNSFP_1=clinvar_id,clinvar_clnsig,clinvar_review,clinvar_trait,1000Gp3_EUR_AF,CADD_phred,Eigen-PC-phred_coding,ExAC_NFE_AF,LRT_pred,
   export dbNSFP_2=FATHMM_pred,GERP++_RS,GTEx_V7_tissue,MutPred_protID,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,SIFT4G_pred,fathmm-MKL_coding_pred,
   export dbNSFP_3=rs_dbSNP151,fathmm-MKL_coding_pred,gnomAD_exomes_NFE_AF,gnomAD_genomes_NFE_AF
   export dbNSFP_fields=${dbNSFP_1}${dbNSFP_2}${dbNSFP_3}
   vep -i ${INF}/work/${s}.vepinput -o ${INF}/work/${s}.dbNSFP --cache --distance 500000 --force --offline --pick --tab \
       --plugin LoF,loftee_path:.,human_ancestor_fa:human_ancestor.fa.gz \
       --plugin dbNSFP,${VEP}/dbNSFP4.0a/dbNSFP4.0a.gz,${dbNSFP_fields}
   cd -
done

cd -

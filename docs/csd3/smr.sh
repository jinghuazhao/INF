#!/usr/bin/bash

# setup
smr --beqtl-summary ~/rds/public_databases/smr/cage_eqtl_data_lite_hg19/CAGE.sparse.lite --extract-snp-p 5e-8 --add-n 2765 --make-besd --out CAGE
smr --beqtl-summary CAGE --show-n 
ln -sf CAGE.besd CAGE_snpid.besd
ln -sf CAGE.epi CAGE_snpid.epi
awk '
{
  chr=$1
  pos=$4
  a1=$5
  a2=$6
  if (a1<a2) snpid="chr" chr ":" pos "_" a1 "_" a2;
  else snpid="chr" chr ":" pos "_" a2 "_" a1;
  $2=snpid
  print
}' OFS='\t' CAGE.esi > CAGE_snpid.esi

# SNP query to INF1.merge.txt
sed '1d' ${INF}/work/INF1.merge | \
cut -f6 | \
sort | \
uniq > INF1.merge.MarkerName
smr --beqtl-summary CAGE_snpid --extract-snp INF1.merge.MarkerName --query 5.0e-8 --out INF1.merge
cut -f5,6 ${INF}/work/INF1.merge | \
awk 'NR>1{print $2,$1}' | \
sort -k1,1 | \
join - <(sed '1d' INF1.merge.txt | sort -k1,1) | \
awk '$2==$11' > INF1.merge-coloc
awk '{print $2,$1}' INF1.merge-coloc | \
sort -k1,2 | \
uniq | \
parallel -j1 -C' ' 'awk -vm={2} -vp={1} "\$5==p && \$6==m" ${INF}/work/INF1.merge'

# formal SMR
echo 6 25392021 33392022 1 > HLA.hg19
sed '1d' ${INF}/work/INF1.merge | \
tr '\t' ' ' | \
parallel -j1 -C' ' '
# PLINK file, INF1.merge.uniprot-prot-markername.gz is generated by merge.sh
  gunzip -c ${INF}/work/INF1.merge.*-{5}-{6}.gz | \
  cut -f3 | \
  sed "1d" > {5}-{6}.snps
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --chr {8} --exclude range HLA.hg19 --extract {5}-{6}.snps --make-bed --out {5}-{6}
# BLD format
  smr --bfile {5}-{6} --make-bld --r --ld-wind 4000 --out {5}-{6}
# SMR and HEIDI test
  smr --bld {5}-{6} --gwas-summary ${INF}/work/{5}.ma --beqtl-summary CAGE_snpid --out {5}-{6} --thread-num 10 
'

(
# collate SMR results prefixed with prot,snpid,gene,cistrans  
  cat *.smr | \
  head -1 | \
  awk -vOFS="\t" "{print \"prot\",\"MarkerName\",\"gene\",\"cistrans\",\$0}"
  sed '1d' ${INF}/work/INF1.merge.cis.vs.trans | \
  cut -d, -f2,5,10,14 | \
  tr ',' ' ' | \
  parallel -j1 -C' ' '
    if [ -f {1}-{2}.smr ]; then
       sed "1d" {1}-{2}.smr | \
       awk -vprot={1} -vsnpid={2} -vgene={3} -vcistrans={4} -vOFS="\t" "{print prot,snpid,gene,cistrans,\$0}"
    fi
  '
) > INF1.merge.smr

(
# add rsid
  export nexp=$(sed '1d' INF1.merge.smr | wc -l)
  awk 'NR==1 {$1="topSNPid toprsid MarkerName rsid";$2="prot";$9="";print}' INF1.merge.smr | \
  awk '{$1=$1};1'
  join -22 ${INF}/work/INTERVAL.rsid <(sed '1d' INF1.merge.smr | sort -k2) | \
  sort -k10 | \
  join -210 ${INF}/work/INTERVAL.rsid - | \
  awk -vnexp=${nexp} '$10!="NA" && $25<=0.05/nexp && $26>0.01 && $26!="NA"'
) > INF1.merge.coloc
awk '$5==$10' INF1.merge.coloc

export nassoc=$(sed '1d' INF1.merge.coloc | wc -l)
export nprot=$(cut -d' ' -f5 INF1.merge.coloc | sort | uniq | wc -l)
export ngene=$(cut -d' ' -f10 INF1.merge.coloc | sort | uniq | wc -l)
echo $nassoc $nexp $nprot $ngene
Rscript -e "with(read.table('INF1.merge.coloc',header=TRUE),table(cistrans))"

# hgnc_symbol <--> probe
awk '!/BDNF/ && NR>1 {
  if($3=="\"Q8NF90\"") $7="\"FGF5\""; else if($3=="\"Q8WWJ7\"") $7="\"CD6\"";
  print
}' FS='\t' OFS='\t' ${INF}/doc/olink.inf.panel.annot.tsv | \
cut -f3,7 | \
sed '1d;s/"//g' | \
sort -k1,1 | \
join -22 - ${INF}/work/inf1.tmp | \
sort -k2,2 | \
join -12 -25 - <(sort -k5,5 CAGE.epi) > INF1.merge.probe

# locus plot
sed '1d' ${INF}/work/INF1.merge | \
tr '\t' ' ' | \
parallel -j1 -C' ' '
  for probe in $(cut -d" " -f3,5 INF1.merge.probe | cut -d" " -f2 | tr "\n" " ")
  do
    export probe=${probe}
    smr --bfile {5}-{6} --gwas-summary ${INF}/work/{5}.ma --beqtl-summary CAGE_snpid \
        --gene-list ~/rds/public_databases/smr/plot/glist_hg19_refseq.txt --probe ${probe} --probe-wind 500 --plot \
        --out smr
    R --no-save -q <<\ \ \ \ END
      HPC_WORK <- Sys.getenv("HPC_WORK")
      probe <- Sys.getenv("probe")
      source(file.path(HPC_WORK,"smr_1.03_src","plot","plot_SMR.r"))
    # Read the data file in R:
      SMRData = ReadSMRData(file.path("plot",paste0("smr.",probe,".txt")))
      png(file.path("plot",paste0(probe,".png")), res=300, units="cm", width=40, height=20)
      SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.01, plotWindow=500, max_anno_probe=16)
    # smr_thresh: genome-wide significance level for the SMR test.
    # heidi_thresh: threshold for the HEIDI test. The default value is 0.05.
    # cis_wind: size of a window centred around the probe to select cis-eQTLs for plot. The default value is 2000Kb.
    # max_anno_probe: maximum number of probe names to be displayed on the figure. The default value is 16.
      dev.off()
    END
  done
'

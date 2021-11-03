#!/usr/bin/bash

function INTERVAL_eQTLGen_SCALLOP()
{
# init -- actually two versions of RNASeq results below gives the same LTBR.lz
  export rnaseq=tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated.csv
  export rnaseq=tensorqtl_allSNPs_MAF0.005_merged_annotated.csv
  grep -w ${rsid} ${rnaseq}
  zgrep ENSG00000256433 ${INF}/work/ensGtp.txt.gz | \
  cut -f2 | \
  zgrep -f - ${INF}/work/ensemblToGeneName.txt.gz
# LocusZoom plot
  read chr start end < st.tmp
  awk -vFS="," -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1 || ($5==chr && $6>=start && $6<=end && index($0,gene)>0)' ${rnaseq} | \
  tr "," "\t" > LTBR.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal LTBR.lz --delim tab title="INTERVAL-LTBR" \
            --markercol variant_id --pvalcol pval --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=INTERVAL --rundir .
  mv INTERVAL_chr${chr}_${bracket}.pdf INTERVAL-LTBR-cis.pdf
# https://www.eqtlgen.org/trans-eqtls.html
# https://www.eqtlgen.org/cis-eqtls.html
{
  export AF=2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz
  export cis=2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  read chr start end < st.tmp
  (
    gunzip -c ${INF}/work/${AF} | awk -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} 'NR==1||($2==chr && $3>=start && $3<=end)' > eQTLGen.AF
    gunzip -c ${cis} | \
    head -1
    gunzip -c ${cis} | \
    awk -vchr=${chr} -vstart=${start} -vend=${end} -vgene=${gene} '$3==chr && $4>=start && $4<=end && index($0,gene)>0' | \
    sort -k3,3n -k4,4n
  ) > eQTLGen.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal eQTLGen.lz --delim tab title="eQTLGen-LTBR" \
            --markercol SNP --pvalcol Pvalue --chr ${chr} --start ${b1} --end ${b2} \
            --no-date --plotonly --prefix=eQTLGen --rundir .
  mv eQTLGen_chr${chr}_${bracket}.pdf eQTLGen-LTBR-cis.pdf
  read chr start end < st.tmp
  (
    gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
    head -1 | \
    awk '{$1=$3 " " $1;$3="";$12="P"};1' | \
    awk '{$1=$1};1'
    gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
    awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} '$1 == chr && $2 >= start && $2 <= end {$12=10^$12;print}' | \
    sort -k3,3 | \
    join -23 <(awk -vchrom=chr${chr} 'index($0,chrom)>0' INTERVAL.rsid) - | \
    sort -k2,2n -k3,3n | \
    cut -d' ' -f1 --complement
  ) | \
  tr ' ' '\t' > TNFB.lz
  rm -f ld_cache.db
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal TNFB.lz --delim tab title="SCALLOP-TNFB" \
            --markercol MarkerName --pvalcol P --chr ${chr} --start ${b1} --end ${b2} --gwas-cat whole-cat_significant-only \
            --no-date --plotonly --prefix=TNFB --rundir .
  mv TNFB_chr${chr}_${bracket}.pdf SCALLOP-TNFB-cis.pdf
}

function stack_assoc_plot_hyprcoloc()
{
  join <(sed '1d' ${INF}/MS/v3.lz | awk '{print $1,$7/$8,$5,$6,$2}' | sort -k1,1) \
       <(sed '1d' ${INF}/work/eQTLGen.lz | \
         awk '{
                 chr=$3;pos=$4;A1=toupper($5);A2=toupper($6);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$7,A1,A2,$2,$3,$4
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10/$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$5,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/work/LTBR.gassoc

  cut -d' ' -f1 ${INF}/work/LTBR.gassoc > ${INF}/work/LTBR.snpid
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/work/LTBR.snpid --r square --out ${INF}/work/LTBR
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/work/LTBR.snpid --freq --out ${INF}/work/LTBR

  Rscript -e '
  INF <- Sys.getenv("INF")
  library(gassocplot)
  d <- read.table(file.path(INF,"work","LTBR.gassoc"),col.names=c("snpid","marker","chr","pos","A1","A2","MS","LTBR","TNFB"))
  markers <- d[c("marker","chr","pos")]
  ld <- read.table(file.path(INF,"work","LTBR.ld"),col.names=with(d,marker),row.names=with(d,marker))
  z <- d[c("MS","LTBR","TNFB")]
  rownames(z) <- with(d,marker)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("MS","LTBR","TNFB"), ylab = "-log10(P)", top.marker="rs1800693",legend=TRUE)
  pdf(file.path(INF,"plots","LTBR.pdf"),width=8,height=13)
  grid::grid.draw(sap)
  dev.off()
# stack_assoc_plot_save(sap, "LTBR.png", 5, width=8, height=13, dpi=300)
  '
  join <(sed '1d' ${INF}/MS/v3.lz | awk '{print $1,$7,$5,$6,$2}' | sort -k1,1) \
       <(Rscript -e '
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim("eQTLGen.lz") %>% filter(GeneSymbol=="LTBR"),read.delim("eQTLGen.AF"),by="SNP") %>%
              mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","b")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$5,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/work/LTBR.beta
  join <(sed '1d' ${INF}/MS/v3.lz | awk '{print $1,$8,$5,$6,$2}' | sort -k1,1) \
       <(Rscript -e '
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim("eQTLGen.lz") %>% filter(GeneSymbol=="LTBR"),read.delim("eQTLGen.AF"),by="SNP") %>%
              mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","se")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    print $1,$5,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/work/LTBR.se
  Rscript -e '
    options(width=200)
    id <- c("marker","chr","pos")
    traits <- c("MS","LTBR","TNFB")
    d <- read.table("LTBR.beta",col.names=c("snpid",id,"A1","A2",traits))
    markers <- d[id]
    betas <- as.matrix(d[traits])
    rownames(betas) <- with(d,marker)
    d <- read.table("LTBR.se",col.names=c("snpid",id,"A1","A2",traits))
    ses <- as.matrix(d[traits])
    rownames(ses) <- with(d,marker)
    hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=with(markers,marker))
  '
}

function PWCoCo()
{
  join <(sed '1d' ${INF}/MS/rs1800693/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$11/$12,$4,$4,$2}' | sort -k1,1) \
       <(sed '1d' ${INF}/work/eQTLGen.lz | \
         awk '{
                 chr=$3;pos=$4;A1=toupper($5);A2=toupper($6);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$7,A1,A2,$2,$3,$4
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10/$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs1800693/LTBR.gassoc

  cut -d' ' -f1 ${INF}/MS/rs1800693/LTBR.gassoc > ${INF}/MS/rs1800693/LTBR.snpid
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/MS/rs1800693/LTBR.snpid --r square --out ${INF}/MS/rs1800693/LTBR
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/MS/rs1800693/LTBR.snpid --freq --out ${INF}/MS/rs1800693/LTBR

  Rscript -e '
  INF <- Sys.getenv("INF")
  library(gassocplot)
  d <- read.table(file.path(INF,"MS","rs1800693","LTBR.gassoc"),col.names=c("snpid","marker","chr","pos","A1","A2","MS","LTBR","TNFB"))
  markers <- d[c("marker","chr","pos")]
  ld <- read.table(file.path(INF,"MS","rs1800693","LTBR.ld"),col.names=with(d,marker),row.names=with(d,marker))
  z <- d[c("MS","LTBR","TNFB")]
  rownames(z) <- with(d,marker)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("MS","LTBR","TNFB"), ylab = "-log10(P)", top.marker="rs2364485",legend=TRUE)
  pdf(file.path(INF,"MS","rs1800693","LTBR.pdf"),width=8,height=13)
  grid::grid.draw(sap)
  dev.off()
  '
  join <(sed '1d' ${INF}/MS/rs1800693/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$11,$4,$4,$2}' | sort -k1,1) \
       <(Rscript -e '
              INF <- Sys.getenv("INF")
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim(file.path(INF,"work","eQTLGen.lz")) %>%
                          filter(GeneSymbol=="LTBR"),read.delim(file.path(INF,"work","eQTLGen.AF")),by="SNP") %>%
                          mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","b")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs1800693/LTBR.beta
  join <(sed '1d' ${INF}/MS/rs1800693/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$12,$4,$4,$2}' | sort -k1,1) \
       <(Rscript -e '
              INF <- Sys.getenv("INF")
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim(file.path(INF,"work","eQTLGen.lz")) %>%
                          filter(GeneSymbol=="LTBR"),read.delim(file.path(INF,"work","eQTLGen.AF")),by="SNP") %>%
                          mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","se")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs1800693/LTBR.se
  Rscript -e '
    options(width=200)
    INF <- Sys.getenv("INF")
    source(file.path(INF,"rsid","LTBR.R"))
    id <- c("marker","chr","pos")
    traits <- c("MS","LTBR","TNFB")
    d <- read.table(file.path(INF,"MS/rs1800693","LTBR.beta"),col.names=c("snpid",id,"A1","A2",traits))
    markers <- d[id]
    print(cor(d[traits]))
    LTBR(d[c("pos","MS","LTBR","TNFB")],file.path(INF,"MS","rs1800693","LTBR.png"))
    betas <- as.matrix(d[traits])
    rownames(betas) <- with(d,marker)
    d <- read.table(file.path(INF,"MS/rs1800693","LTBR.se"),col.names=c("snpid",id,"A1","A2",traits))
    ses <- as.matrix(d[traits])
    rownames(ses) <- with(d,marker)
    hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=with(markers,marker))
  '
# --- rs2364485
  join <(sed '1d' ${INF}/MS/rs2364485/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$11/$12,$4,$4,$2}' | sort -k1,1) \
       <(sed '1d' ${INF}/work/eQTLGen.lz | \
         awk '{
                 chr=$3;pos=$4;A1=toupper($5);A2=toupper($6);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$7,A1,A2,$2,$3,$4
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10/$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs2364485/LTBR.gassoc

  cut -d' ' -f1 ${INF}/MS/rs2364485/LTBR.gassoc > ${INF}/MS/rs2364485/LTBR.snpid
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/MS/rs2364485/LTBR.snpid --r square --out ${INF}/MS/rs2364485/LTBR
  plink --bfile ${INF}/INTERVAL/cardio/INTERVAL --extract ${INF}/MS/rs2364485/LTBR.snpid --freq --out ${INF}/MS/rs2364485/LTBR
  Rscript -e '
  INF <- Sys.getenv("INF")
  library(gassocplot)
  d <- read.table(file.path(INF,"MS","rs2364485","LTBR.gassoc"),col.names=c("snpid","marker","chr","pos","A1","A2","MS","LTBR","TNFB"))
  markers <- d[c("marker","chr","pos")]
  ld <- read.table(file.path(INF,"MS","rs2364485","LTBR.ld"),col.names=with(d,marker),row.names=with(d,marker))
  z <- d[c("MS","LTBR","TNFB")]
  rownames(z) <- with(d,marker)
  sap <- stack_assoc_plot(markers, z, ld, traits = c("MS","LTBR","TNFB"), ylab = "-log10(P)", top.marker="rs1800693",legend=TRUE)
  pdf(file.path(INF,"MS","rs2364485","LTBR.pdf"),width=8,height=13)
  grid::grid.draw(sap)
  dev.off()
  '
  join <(sed '1d' ${INF}/MS/rs2364485/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$11,$4,$4,$2}' | sort -k1,1) \
       <(Rscript -e '
              INF <- Sys.getenv("INF")
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim(file.path(INF,"work","eQTLGen.lz")) %>%
                          filter(GeneSymbol=="LTBR"),read.delim(file.path(INF,"work","eQTLGen.AF")),by="SNP") %>%
                          mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","b")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$10,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    if($3!=$7) $6=-$6
    if($3!=$13) $12=-$12
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs2364485/LTBR.beta
  join <(sed '1d' ${INF}/MS/rs2364485/EUR-v3.cma.cojo | awk '{if($12!="NA") print $2,$12,$4,$4,$2}' | sort -k1,1) \
       <(Rscript -e '
              INF <- Sys.getenv("INF")
              suppressMessages(library(dplyr))
              cis_pQTL <- merge(read.delim(file.path(INF,"work","eQTLGen.lz")) %>%
                          filter(GeneSymbol=="LTBR"),read.delim(file.path(INF,"work","eQTLGen.AF")),by="SNP") %>%
                          mutate(data.frame(gap::get_b_se(AlleleB_all,NrSamples,Zscore)))
              write.table(cis_pQTL[c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","se")], sep = "\t",
                          row.names = FALSE, col.names = TRUE, quote=FALSE)
             ' | \
         awk '{
                 chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                 if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                 print snpid,$6,A1,A2,$1,$2,$3
              }' | \
         sort -k1,1 \
         ) | \
  join - <(sed '1d' ${INF}/work/TNFB.lz | \
           awk '{
                   chr=$2;pos=$3;A1=toupper($4);A2=toupper($5);
                   if(A1<A2) snpid="chr"chr":"pos"_"A1"_"A2;else snpid="chr"chr":"pos"_"A2"_"A1;
                   print snpid,$11,A1,A2,$1
                }' | \
           sort -k1,1 \
          ) | \
  awk -vOFS="\t" '
  {
    print $1,$9,$10,$11,$3,$4,$2,$6,$12
  }' | \
  awk 'a[$1]++==0' | \
  awk 'a[$2]++==0' | awk '$2!="NA"' > ${INF}/MS/rs2364485/LTBR.se
  Rscript -e '
    options(width=200)
    INF <- Sys.getenv("INF")
    source(file.path(INF,"rsid","LTBR.R"))
    id <- c("marker","chr","pos")
    traits <- c("MS","LTBR","TNFB")
    d <- read.table(file.path(INF,"MS/rs2364485","LTBR.beta"),col.names=c("snpid",id,"A1","A2",traits))
    markers <- d[id]
    print(cor(d[traits]))
    LTBR(d[c("pos","MS","LTBR","TNFB")],file.path(INF,"MS","rs2364485","LTBR.png"))
    betas <- as.matrix(d[traits])
    rownames(betas) <- with(d,marker)
    d <- read.table(file.path(INF,"MS/rs2364485","LTBR.se"),col.names=c("snpid",id,"A1","A2",traits))
    ses <- as.matrix(d[traits])
    rownames(ses) <- with(d,marker)
    hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=with(markers,marker))
  '
}

cd work
export chr=12
export pos=6514963
export gene=LTBR
export rsid=rs2364485
export flank_kb=1000
export b1=6300000
export b2=6700000
export bracket=${b1}-${b2}

tabix ${INF}/METAL/gwas2vcf/TNFB.tsv.gz ${chr}:${bracket} > TNFB.tbx

module load python/2.7
awk -vchr=${chr} -vpos=${pos} -vd=$((${flank_kb}*1000)) 'BEGIN{print chr,pos-d,pos+d}' > st.tmp

INTERVAL_eQTLGen_SCALLOP
stack_assoc_plot_hyprcoloc
PWCoCo
cd -

# MS v3/v3.lz
# CHR BP SNP A1 A2 N P OR
# SNPid   SNP     chr     pos     a1      a2      b       se      p
# eQTLGen.lz
# Pvalue	SNP	SNPChr	SNPPos	AssessedAllele	OtherAllele	Zscore	Gene	GeneSymbol	GeneChr	GenePos	NrCohorts	NrSamples	FDR	BonferroniP
# TNFB.lz
# MarkerName	Chromosome	Position	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Effect	StdErr	P	Direction	HetISq	HetChiSq	HetDf	logHetP	N

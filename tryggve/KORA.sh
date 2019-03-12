#!/usr/bin/bash
# 12-3-2019 JHZ

module load bcftools/1.9
module load plink2/1.90beta5.4
module load bolt-lmm/2.3.2
module load gcc/5.4.0 lapack/3.8.0 qctool/2.0.1
module load intel/redist/2019 intel/perflibs/64/2019 gcc/5.4.0 R/3.5.0-ICC-MKL
module load snptest/2.5.2

# module load gcta/1.91.0beta
# but it doesn't handle --grm and --out share the same file name
# so /data/jinhua/gcta_1.91.7beta/gcta64 is called through symbolic link

export rt=$HOME/INF
cd $rt/KORA

function info()
{
  ## concatenate .info and filter .info>=0.4
  export info=/data/jinhua/data/KORA/impute_out_info/info

  ls $info | \
  sed 's/KORAS4F4_AffyAxiom_N3788_chr//g;s/_imputation.impute2_info\*//g;s/_/ /g' | \
  sort -t' ' -k1,1g -k2,2g | \
  parallel -j1 --env info -C' '  '
  (
    head -1 $info/KORAS4F4_AffyAxiom_N3788_chr{1}_1_imputation.impute2_info
    awk "NR>1" $info/KORAS4F4_AffyAxiom_N3788_chr{1}_{2}_imputation.impute2_info
  ) > chr{1}.info
  '
  seq 22 | \
  parallel -j1 -C' ' 'awk -vchr={} "NR>1 && \$7>=0.4 {print chr \":\" \$3}" chr{}.info > chr{}.id
# only ~1M variants left with info>=0.4, so preferably skipped
  plink --vcf chr{}.vcf.gz --extract chr{}.id --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out chr{}
  plink --vcf chr{}.vcf.gz --extract chr{}.prune.in --make-bed --out chr{}'
}

function ln()
{
  ## SNP filtering for .vcf.gz
  export vcf="/data/jinhua/data/KORA/Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220"

  seq 22 | \
  parallel -j3 -C' ' 'ln -sf Affy\ AxiomPhase3_n3775_CodeAX1KG3_V2_LU9220/chr{}_N3775_1000GPhase3.dose.vcf.gz chr{}.vcf.gz'
}

function id()
{
  bcftools query -l chr22.vcf.gz | \
  sort > genotype.id
  (
    echo ID
    echo 0
    cat genotype.id
  ) > KORA.samples
  export llod=/data/jampet/KORA/kora.below.llod.normalised.prot.txt
  export llod=/data/jampet/KORA/kora.below.llod.normalised.prot.2019-03-08.txt
  awk -vOFS="\t" '{
    $1=$1 OFS $1
    if($3=="M") $3=0; else if($3=="F") $3=1
    print
  }' $llod | \
  awk -vOFS="\t" '{if(NR==1) {$1="FID"; $2="IID"}};1' > phenocovar.txt
  cut -f1-2 phenocovar.txt | \
  awk 'NR>1' | \
  sort -k1,1 | \
  join genotype.id - | \
  cut -d' ' -f1 > protein.id
  join -v1 genotype.id protein.id | \
  awk -vOFS="\t" '{print $1,$1}' > remove.id
}

function snp()
{
  echo "--> -sample-stats and -snp-stats"
  seq 22 | \
  parallel -j12 'bcftools view -S protein.id -O z chr{}.vcf.gz -o KORA{}.vcf.gz'
  qctool -g KORA#.vcf.gz -s KORA.samples -excl-samples remove.id -vcf-genotype-field GP \
         -sample-stats -osample KORA.sample-stats -snp-stats -osnp KORA#.snp-stats -threads 12
  echo "--> PCs based on independent SNPs"
  seq 22 | \
  parallel -j1 -C' ' '
    bcftools annotate --set-id "chr%CHROM\:%POS\_%REF\_%ALT" chr{}.vcf.gz -O z -o KORA{}.vcf.gz
    plink --vcf KORA{}.vcf.gz --list-duplicate-vars --out chr{}
    awk "NR>1{split(\$NF,dupids,\" \");print dupids[1]}" chr{}.dupvar > chr{}.dupid
    bcftools query -i "MAF>0.01 && R2>=0.4" -f"%ID\n" KORA{}.vcf.gz | \
    join -v2 chr{}.dupid - > chr{}.mafr2
    plink --vcf KORA{}.vcf.gz --extract chr{}.mafr2 --remove remove.id --make-bed --out nodup{}
    awk -vOFS="\t" "
    {
      CHR=\$1
      POS=\$4
      a1=\$5
      a2=\$6
      if (a1>a2) snpid=\"chr\" CHR \":\" POS \"_\" a2 \"_\" a1;
      else snpid=\"chr\" CHR \":\" POS \"_\" a1 \"_\" a2
      print snpid, \$2
    }" nodup{}.bim > nodup{}.snpid
    plink --bfile nodup{} --update-name nodup{}.snpid 1 2 --make-bed --out KORA{}
  '
  seq 22 | \
  awk -vp=KORA '{print p NR}' > KORA.list
  plink --merge-list KORA.list --make-bed --out KORA
  plink --bfile KORA --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out KORA
  plink --bfile KORA --extract KORA.prune.in --make-bed --out KORA.prune
  gcta64 --bfile KORA.prune --make-grm-bin --thread-num 5 --out KORA
  gcta64 --grm KORA --pca 5 --out KORA
  echo "--> relatedness"
  king -b KORA.prune.bed --related --prefix KORA.prune
  awk 'NR>1{print $4}' KORA.prune.kin0 > KORA.prune.relatedness
}

function phenocovar()
{
  R -q --no-save <<\ \ END
    phenocovar <- read.delim("phenocovar.txt",as.is=TRUE)
    sample_stats <- read.delim("KORA.sample-stats",skip=12,nrows=1070,as.is=TRUE)
    missing_proportion <- with(sample_stats,{data.frame(FID=sample,IID=sample,missing=missing_proportion)})
    eigenvec <- read.table("KORA.eigenvec",col.names=c("FID","IID",paste0("PC",1:5)))
    PCs <- paste0("PC",1:5)
    eigenvec[PCs] <- eigenvec[PCs]*100
    pheno <- merge(missing_proportion,merge(eigenvec,phenocovar,by=c("FID","IID")),by=c("FID","IID"))
    names(pheno)[1:2] <- c("ID_1","ID_2")
    l2 <- c(rep("0",3),rep("C",5+2),rep("P",91))
    write.table(rbind(l2,pheno),file="KORA.pheno",quote=FALSE,row.names=FALSE)
  END
  cut -f5-95 phenocovar.txt | \
  awk 'NR==1{gsub(/UH_O_/,"");gsub(/\t/," ");print}' > KORA.varlist
}

function snptest_assoc()
{
  parallel -j12 --env rt -C' ' '
    snptest \
      -data KORA{2}.vcf.gz KORA.pheno \
      -exclude_samples KORA.prune.relatedness \
      -o {1}-{2} \
      -printids \
      -lower_sample_limit 50 \
      -frequentist 1 \
      -genotype_field GP \
      -missing_code NA,-999 \
      -method expected \
      -pheno UH_O_{1} \
      -cov_all \
      -use_raw_covariates \
      -use_raw_phenotypes \
      -use_long_column_naming_scheme \
      -hwe \
      -log {1}-{2}.log;gzip -f {1}-{2};gzip -f {1}-{2}.log' ::: $(cat KORA.varlist) ::: $(seq 22)
  for p in $(cat KORA.varlist)
  do
    echo $p
    export prot=$p
    seq 22 | \
    parallel -j1 --env prot -C' ' '
    (
      gunzip -c ${prot}-{}.gz | \
      awk "NR>21"
    )' | \
    grep -v -E 'not|Completed' | \
    awk 'NR==1 || $3!="chromosome"' | \
    gzip -f > ${p}.gz
  done
}

function bolt_assoc()
## association analysis with BOLT but abandoned for many failed runs
# https://data.broadinstitute.org/alkesgroup/BOLT-LMM/#x1-220005.1.2
{
  seq 22 | \
  awk '{print "protein" $1 ".gen.gz"}' > bolt.list
  awk -vOFS="\t" '{print $1, $1}' protein.id > bolt.id
  echo "--> IMPUTE2 format"
  seq 22 | \
  parallel -j5 -C' ' 'bcftools convert --samples-file protein.id --tag GP chr{}.vcf.gz -g protein{}'
  seq 22 | \
  parallel -j2 -C' ' '
  bolt \
    --bfile KORA.prune \
    --impute2FileList=bolt.list \
    --impute2FidIidFile=bolt.id \
    --LDscoresUseChip \
    --maxModelSnps 50000000 \
    --noMapCheck \
    --phenoFile=phenocovar.txt --phenoCol UH_O_{} \
    --covarFile=phenocovar.txt --covarCol sex --qCovarCol age \
    --remove remove.id \
    --lmm --statsFile={}.stats 2>&1 | \
  tee {}.log' ::: $(cat KORA.varlist)
}

function qqman()
{
  cat $rt/sumstats/KORA.list | \
  parallel -j5 -C' ' '
      export gene={1}
      export protein={3}
      R --no-save -q <<\ \ \ \ \ \ END
      gene <- Sys.getenv("gene");
      print(gene);
      protein <- Sys.getenv("protein");
      print(protein);
      gz <- gzfile(paste0(gene,".gz"));
      .libPaths("/services/tools/R/3.5.0/lib64/R/library")
      require(qqman);
      tbl <- read.table(gz,as.is=TRUE,header=TRUE);
      tbl <- within(tbl,{
         SNP <- rsid 
         CHR <- as.numeric(unlist(strsplit(rsid,":"))[1])
         BP <- position
         P <- tbl[[paste0("UH_O_",gene,"_PC1_PC2_PC3_PC4_PC5_age_sex_frequentist_add_expected_pvalue")]]
      })
      tbl <- subset(tbl,!is.na(CHR)&!is.na(BP)&!is.na(P))
      qq <- paste0(protein,".qq.png");
      png(qq,width=12,height=10,units="in",pointsize=4,res=300)
      qq(with(tbl,P))
      dev.off()
      manhattan <- paste0(protein,".manhattan.png");
      png(manhattan,width=12,height=10,units="in",pointsize=4,res=300)
      manhattan(tbl,main=protein,genomewideline=-log10(5e-10),suggestiveline=FALSE,ylim=c(0,25));
      dev.off();
      END
  '
}

function h2()
{
  awk 'NR>1{$3="";$4="";print}' phenocovar.txt | \
  awk '{$1=$1;print}' > prot.dat
  awk 'NR>1{print $1,$2,$3}' phenocovar.txt > age.dat
  awk 'NR>1{print $1,$2,$4}' phenocovar.txt > sex.dat
  awk 'NR==1{gsub(/\t/, "\n", $0); print}' phenocovar.txt | \
  awk 'NR>4' | \
  awk '{gsub(/UH_O_/,"",$1);print $1,NR}' | \
  parallel -j1 -C' ' '
    gcta64 --reml --grm KORA --pheno prot.dat --mpheno {2} --covar sex.dat --qcovar age.dat \
           --thread-num 5 --out {1} 2>&1 | \
    tee {1}.log
  '
  grep '/' *hsq | \
  sed 's|.hsq:V(G)/Vp||g' > h2.out
  (
    echo "prot h2 se p"
    grep 'Pval' *hsq | \
    sed 's/.hsq:Pval//g' | \
    join h2.out - 
  ) > h2.stats
}

snptest_assoc

#/usr/bin/bash

module load gcc/6
if [ ! -d ${INF}/mr/pQTLs ]; then mkdir -p ${INF}/mr/pQTLs; fi

function iv()
{
  (
    echo SNP Phenotype effect_allele other_allele eaf beta se pval N
  # rsid prot Allele1 Allele2 Freq1 Effect StdErr log.P. cis.trans
    cut -f2,3,6,7,8-11,17,21 ${INF}/work/INF1.METAL | \
    awk -vtype=${type} '$10==type {print $1,$2,toupper($3),toupper($4),$5,$6,$7,10^$8,$9}'
  ) > ${INF}/mr/pQTLs/INF1_${type}.ins
}

function INF1_efo()
{
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    outcomes <- c("ieu-a-7","ebi-a-GCST007432")
    ieugwasr::gwasinfo(id = outcomes)
    type <- Sys.getenv("type")
    ivs <- read.table(file.path(INF,"mr","pQTLs",paste0("INF1_",type,".ins")),as.is=TRUE,header=TRUE)
    Ins <- TwoSampleMR::format_data(ivs,snp_col="SNP",samplesize_col="N")
    Ids <- TwoSampleMR::extract_outcome_data(snps=with(Ins,SNP),outcomes=outcomes)
    harmonise <- TwoSampleMR::harmonise_data(Ins, Ids)
    prefix <- file.path(INF,"mr","pQTLs",paste0("INF1_pQTL-combined-",type,"-"))
    pQTLtools::pqtlMR(pqtlMRinput=list(Ins,Ids,harmonise),prefix=prefix)
    prefix <- file.path(INF,"mr","pQTLs",paste0("INF1_rev_pQTL-combined-",type,"-"))
    pQTLtools::pqtlMR(pqtlMRinput=list(Ins,Ids,harmonise),prefix=prefix,reverse=TRUE)
  END
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    outcomes <- with(read.delim(file.path(INF,"rsid","efo.txt")),MRBASEID)
    type <- Sys.getenv("type")
    ivs <- read.table(file.path(INF,"mr","pQTLs",paste0("INF1_",type,".ins")),as.is=TRUE,header=TRUE)
    Ins <- TwoSampleMR::format_data(ivs,snp_col="SNP",samplesize_col="N")
    Ids <- TwoSampleMR::extract_outcome_data(snps=with(Ins,SNP),outcomes=outcomes)
    harmonise <- TwoSampleMR::harmonise_data(Ins,Ids)
    prefix <- file.path(INF,"mr","pQTLs",paste0("efo_pQTL-combined-",type,"-"))
    pQTLtools::pqtlMR(pqtlMRinput=list(Ins, Ids, harmonise),prefix=prefix)
    prefix <- file.path(INF,"mr","pQTLs",paste0("efo_rev_pQTL-combined-",type,"-"))
    pQTLtools::pqtlMR(pqtlMRinput=list(Ins, Ids, harmonise),prefix=prefix,reverse=TRUE)
  END
}

function collect()
{
  echo ${prefix} -- ${id} -- ${trait}
  (
    cat ${prefix}result.txt | head -1
    grep -w ${id} ${prefix}result.txt | grep "Wald ratio"
  ) | grep -v _rev_ > ${prefix}${id}.result
  (
    cat ${prefix}single.txt | head -1
    grep -w ${id} ${prefix}single.txt | grep -v -e Egger -e Inverse
  ) | grep -v _rev_ > ${prefix}${id}.single
}

function collect_rev()
{
  echo ${prefix} -- ${id} -- ${trait}
  (
    cat ${prefix}result.txt | head -1
    grep -w ${id} ${prefix}result.txt | grep "Wald ratio"
  ) | awk -v FS="\t" -v id=${id} 'NR==1||$1==id' > ${prefix}${id}.result
  (
    cat ${prefix}single.txt | head -1
    grep -w ${id} ${prefix}single.txt | grep -v -e Egger -e Inverse
  ) | awk -v FS="\t" -v id=${id} 'NR==1||$3==id' > ${prefix}${id}.single
}

function collect_all()
{
  for i in ieu-a-7 ebi-a-GCST007432
  do
      export id=${i}
      if [ "${i}" == "ieu-a-7" ]; then
         export trait="CHD || ${i}"
      else
         export trait="FEV1 || ${i}"
      fi
      export prefix=${INF}/mr/pQTLs/INF1_pQTL-combined-${type}-
      collect
      export prefix=${INF}/mr/pQTLs/INF1_rev_pQTL-combined-${type}-
      collect_rev
  done
  export nrows=$(sed '1d' ${INF}/rsid/efo.txt | wc -l | cut -d' ' -f1)
  for i in $(seq ${nrows})
  do
      export trait=$(sed '1d' ${INF}/rsid/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $2}')
      export id=$(sed '1d' ${INF}/rsid/efo.txt | awk -vFS="\t" -vnr=${i} 'NR==nr{print $4}')
      export prefix=${INF}/mr/pQTLs/efo_pQTL-combined-${type}-
      collect
      export prefix=${INF}/mr/pQTLs/efo_rev_pQTL-combined-${type}-
      collect_rev
  done
}

for type in cis trans
do
  export type=${type}
  iv
  INF1_efo
  collect_all
done

# Bidirectionality test for FGF.5
function dummy()
{
R --no-save -q <<END
  options(width=200)
  INF1_cis <- read.delim("INF1_cis.ins",sep=" ")
  pQTLtools::pqtlMR(subset(INF1_cis,Phenotype=="FGF.5"),"ieu-a-7",prefix="test")
  h <- read.delim("test-harmonise.txt")
  require(TwoSampleMR)
  info <- ieugwasr::gwasinfo("ieu-a-7")
# binary outcomes
  lor <- with(h,beta.outcome)
  af <- with(h,eaf.outcome)
  ncase <- with(info,ncase)
  ncontrol <- with(info,ncontrol)
  prevalence <- 0.1
  pval.exposure <- with(h,pval.exposure)
  samplesize.exposure <- 11787
  outcome <- with(info,sample_size)
  r.exposure <- get_r_from_pn(pval.exposure,samplesize.exposure)
  r.outcome <- get_r_from_lor(lor, af, ncase, ncontrol, prevalence, model = "logit", correction = FALSE)
  h <- data.frame(h,samplesize.exposure=11787,r.exposure=r.exposure,r.outcome=r.outcome)
  directionality_test(h)
END
}

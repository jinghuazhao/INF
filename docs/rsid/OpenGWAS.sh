#!/usr/bin/bash

if [ ! -d ${INF}/OpenGWAS ]; then mkdir ${INF}/OpenGWAS; fi
cd ${INF}/OpenGWAS
export CADFEV1=(ebi-a-GCST005195 ebi-a-GCST007432) # CAD does not contain AF
export CADFEV1=(ebi-a-GCST004787 ebi-a-GCST007432) # CAD is 2017 but unavailable!
export CADFEV1=(ebi-a-GCST003116 ebi-a-GCST007432) # CAD is 2015
for id in $(echo $(cut -f4 ${INF}/rsid/efo.txt | sed '1d') ${CADFEV1[*]})
do
  export f=https://gwas.mrcieu.ac.uk/files/${id}/${id}.vcf.gz
  if [ ! -f ${id}.vcf.gz ]; then wget ${f}; fi
  if [ ! -f ${id}.vcf.gz.tbi ]; then wget ${f}.tbi; fi
done
cd -

R --no-save -q <<END
  options(width=200)
  library(dplyr)
  INF <- Sys.getenv("INF")
  HPC_WORK <- Sys.getenv("HPC_WORK")
  efo <- read.delim(file.path(INF,"rsid","efo.txt"))
  cadfev1 <- c("ebi-a-GCST003116", "ebi-a-GCST007432")
  opengwas_ids <- efo[["MRBASEID"]] %>%
                  c(cadfev1)
  opengwas_info <- ieugwasr::gwasinfo(opengwas_ids)
  write.table(opengwas_info,file=file.path(INF,"OpenGWAS","ieu.info"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(opengwas_info%>%select(id, sample_size, ncontrol, ncase),
              file=file.path(INF,"OpenGWAS","ieu.N"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
  unavail <- c("ieu-b-18","finn-a-M13_SLE","finn-a-D3_SARCOIDOSIS")
  opengwas_ids <- subset(opengwas_ids,!opengwas_ids%in%unavail)
  region <- "1:100-2000000"
  gwasvcf::set_bcftools(path=file.path(HPC_WORK,"bin","bcftools"))
  library(pQTLtools)
  summary_list = purrr::map(opengwas_ids[1:2], ~import_OpenGWAS(., region))
END

# ieu-b-18.vcf.gz
# finn-a-M13_SLE.vcf.gz
# finn-a-D3_SARCOIDOSIS.vcf.gz
# CAD/FEV1, ebi-a-GCST003116

export EBI=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004787/
export GCST004787=${EBI}/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz
wget ${GCST004787} -O ${INF}/OpenGWAS/GCST004787.gz
# It does not have case/controal classification
export CAD2017=~/rds/results/public/gwas/chd/Nelson_2017/raw_results/ukbiobank/ukbb.1000g.exome.soft.meta.filter.out.gz

R --no-save -q <<END
  get_studies(study_id = 'GCST004787')
  get_associations(study_id = 'GCST004787')
END

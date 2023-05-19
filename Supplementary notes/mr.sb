#!/usr/bin/bash

#SBATCH --job-name=gsmr
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --array=1-58
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/mr/gsmr/slurm/_gsmr_%A_%a.out
#SBATCH --error=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/mr/gsmr/slurm/_gsmr_%A_%a.err
#SBATCH --export ALL

export TMPDIR=/rds/user/jhz22/hpc-work/work
export job=$SLURM_ARRAY_TASK_ID

. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4

function rmdup()
# Remove duplicates via PLINK
{
  echo ${chr} ${start} ${end} > ${INF}/mr/gsmr/ref/${prot}-cis.bed1
  plink2 --pfile ${pgen}/impute_dedup_${chr}_interval --extract bed1 ${INF}/mr/gsmr/ref/${prot}-cis.bed1 \
         --make-bed --rm-dup force-first list --out ${INF}/mr/gsmr/ref/${prot}-cis
}

function rmall()
# Remove all duplicates by brute force
{
    awk -vchr=${chr} -vstart=${start} -vend=${end} '$1==chr && $2>=start && $2<end {print $3}' ${pgen}/impute_dedup_${chr}_interval.pvar \
        > ${INF}/mr/gsmr/ref/${prot}-cis.list
    join -v1 <(sort ${INF}/mr/gsmr/ref/${prot}-cis.list) \
             <(awk 'a[$1]++>0' ${INF}/mr/gsmr/ref/${prot}-cis.list | sort) \
             > ${INF}/mr/gsmr/ref/${prot}-cis.keep
    plink2 --pfile ${pgen}/impute_dedup_${chr}_interval \
           --extract ${INF}/mr/gsmr/ref/${prot}-cis.keep \
           --make-bed --out ${INF}/mr/gsmr/ref/${prot}-cis
}

function bfile()
# reference with rsids
{
  if [ ! -d ${INF}/mr/gsmr/ref ]; then mkdir ${INF}/mr/gsmr/ref; fi
  export pgen=~/rds/post_qc_data/interval/imputed/uk10k_1000g_b37/imputed/plink_format/pgen
  awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  awk -vM=1e6 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2, $3, $4-M, $5+M}' | \
  while read prot chr start end
  do
    export prot=${prot}
    export chr=${chr}
    export start=${start}
    export end=${end}
    rmdup
  done
}

function gsmr()
# first list of traits
{
  if [ ! -d ${INF}/mr/gsmr/out ]; then mkdir -p ${INF}/mr/gsmr/out; fi
  awk -vFS="\t" 'NR>1 {print $4}' ${INF}/rsid/efo.txt | \
  while read efo
  do
    export trait=${efo}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' | \
    while read prot
    do
      export prot=${prot}
      echo ${INF}/mr/gsmr/ref/${prot} > ${INF}/mr/gsmr/ref/gsmr_${prot}
      echo ${prot} ${INF}/mr/gsmr/prot/${prot}.gz > ${INF}/mr/gsmr/prot/gsmr_${prot}
      echo ${trait} ${INF}/mr/gsmr/trait/${trait}-${prot}.gz > ${INF}/mr/gsmr/trait/gsmr_${trait}-${prot}

      gcta-1.9 --mbfile ${INF}/mr/gsmr/ref/gsmr_${prot} \
               --gsmr-file ${INF}/mr/gsmr/prot/gsmr_${prot} ${INF}/mr/gsmr/trait/gsmr_${trait}-${prot} \
               --gsmr-direction 0 \
               --clump-r2 0.1 --gwas-thresh 5e-8 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 10 --effect-plot \
               --out ${INF}/mr/gsmr/out/${trait}-${prot}
      if [ ! -f ${INF}/mr/gsmr/out/${trait}-${prot}.eff_plot.gz ]; then continue; fi
      R --no-save -q <<\ \ \ \ \ \ END
        INF <- Sys.getenv("INF")
        trait <- Sys.getenv("trait")
        p <- Sys.getenv("prot")
        source(file.path(INF,"rsid","gsmr_plot.r"))
        eff_dat <- paste0(INF,"/mr/gsmr/out/",trait,"-",p,".eff_plot.gz")
        gsmr_data <- read_gsmr_data(eff_dat)
        gsmr_summary(gsmr_data)
        pdf(paste0(INF,"/mr/gsmr/out/",trait,"-",p,".eff_plot.pdf"))
        par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
        plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
        dev.off()
      END
    done
  done
}

function gsmr2()
# revised list without ukb, bbj, fenn
{
  if [ ! -d ${INF}/mr/gsmr/out ]; then mkdir -p ${INF}/mr/gsmr/out; fi
  cat <(awk -vFS="\t" 'NR>1 && $4 !~ /ukb/ {print $4}' ${INF}/rsid/efo.txt) \
      <(awk -vFS="\t" '/ebi|ieu|finn/{print $4}' ${INF}/OpenGWAS/ukb-replacement.txt) | \
  while read efo
  do
    export trait=${efo}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' | \
    while read prot
    do
      export prot=${prot}
      echo ${INF}/mr/gsmr/ref/${prot} > ${INF}/mr/gsmr/ref/gsmr_${prot}
      echo ${prot} ${INF}/mr/gsmr/prot/${prot}.gz > ${INF}/mr/gsmr/prot/gsmr_${prot}
      echo ${trait} ${INF}/mr/gsmr/trait/${trait}-${prot}.gz > ${INF}/mr/gsmr/trait/gsmr_${trait}-${prot}

      gcta-1.9 --mbfile ${INF}/mr/gsmr/ref/gsmr_${prot} \
               --gsmr-file ${INF}/mr/gsmr/prot/gsmr_${prot} ${INF}/mr/gsmr/trait/gsmr_${trait}-${prot} \
               --gsmr-direction 0 \
               --clump-r2 0.1 --gwas-thresh 5e-8 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 10 --effect-plot \
               --out ${INF}/mr/gsmr/out/${trait}-${prot}
      if [ ! -f ${INF}/mr/gsmr/out/${trait}-${prot}.eff_plot.gz ]; then continue; fi
      R --no-save -q <<\ \ \ \ \ \ END
        INF <- Sys.getenv("INF")
        trait <- Sys.getenv("trait")
        p <- Sys.getenv("prot")
        source(file.path(INF,"rsid","gsmr_plot.r"))
        eff_dat <- paste0(INF,"/mr/gsmr/out/",trait,"-",p,".eff_plot.gz")
        gsmr_data <- read_gsmr_data(eff_dat)
        gsmr_summary(gsmr_data)
        pdf(paste0(INF,"/mr/gsmr/out/",trait,"-",p,".eff_plot.pdf"))
        par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
        plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
        dev.off()
      END
    done
  done
}

function gsmr3()
# revised list without ukb, bbj, fenn AND trait based on extract_outcome_data but now VCF
{
  if [ ! -d ${INF}/mr/gsmr/out ]; then mkdir -p ${INF}/mr/gsmr/out; fi
  export EFO_UPDATE=${INF}/OpenGWAS/efo-update.txt
# cat <(awk -vFS="\t" 'NR>1 && $4 !~ /ukb/ {print $4}' ${INF}/rsid/efo.txt) \
#     <(awk -vFS="\t" '/ebi|ieu|finn/{print $4}' ${INF}/OpenGWAS/ukb-replacement.txt) | \
# cat <(sed '1d' ${EFO_UPDATE} | grep -v -e finn -e ebi-a-GCST90014325 | awk -vFS="\t" '{print $6,2/(1/$2+1/$3)}') \
#     <(sed '1d' ${EFO_UPDATE} | grep -e finn -e ebi-a-GCST90014325 | awk -vFS="\t" '{print $6,2/(1/$2+1/$3)}')  | \
  cat <(sed '1d' ${EFO_UPDATE} | awk -vFS="\t" '{print $1,2/(1/$3+1/$4)}') | \
  while IFS=$'\t' read -r efo trait ncase ncontrol
  do
    export efo=${efo}
    export trait=${trait}
    awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
    awk -vM=1e6 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}' | \
    while read prot
    do
      export prot=${prot}
      export rsid=$(awk -vprot=${prot} '$3==prot && $21=="cis" {print $2}' ${INF}/work/INF1.METAL)
      echo ${INF}/mr/gsmr/ref/${prot} > ${INF}/mr/gsmr/ref/gsmr_${prot}
      echo ${prot} ${INF}/mr/gsmr/prot/${prot}.gz > ${INF}/mr/gsmr/prot/gsmr_${prot}
    # echo ${trait} ${INF}/mr/gsmr/trait/${prot}-${trait}.gz > ${INF}/mr/gsmr/trait/gsmr_${prot}-${trait}
      grep -w ${prot} ${INF}/TNFB/cis.dat | cut -d' ' -f4-6 | while read chr start end
      do
        module load python/2.7
        mkdir ${INF}/mr/gsmr/trait/${prot}-${efo}
        locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${INF}/mr/gsmr/trait/${prot}-${efo}-rsid.txt \
                  --delim space title="${prot}-${efo}-${trait}" \
                  --markercol SNP --pvalcol p --chr ${chr} --start ${start} --end ${end} --cache None \
                  --no-date --plotonly --prefix=${prot}-${efo} --rundir ${INF}/mr/gsmr/trait/${prot}-${efo} --svg --refsnp ${rsid}
      done
      gcta-1.9 --mbfile ${INF}/mr/gsmr/ref/gsmr_${prot} \
               --gsmr-file ${INF}/mr/gsmr/prot/gsmr_${prot} ${INF}/mr/gsmr/trait/gsmr_${prot}-${efo} \
               --gsmr-direction 0 \
               --clump-r2 0.1 --gwas-thresh 5e-8 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 3 --effect-plot \
               --out ${INF}/mr/gsmr/out/${prot}-${efo}
      if [ ! -f ${INF}/mr/gsmr/out/${prot}-${efo}.eff_plot.gz ]; then continue; fi
      gunzip -c ${INF}/mr/gsmr/out/${prot}-${efo}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect | cut -d' ' -f1 | \
      grep -f - ${INF}/work/INTERVAL.rsid > ${INF}/mr/gsmr/out/${prot}-${efo}.rsid
      R --no-save -q <<\ \ \ \ \ \ END
        INF <- Sys.getenv("INF")
        efo <- Sys.getenv("efo")
        trait <- Sys.getenv("trait")
        p <- Sys.getenv("prot")
        r <- Sys.getenv("rsid")
        rt <- paste0(INF,"/mr/gsmr/out/",p,"-",efo)
        snpid_rsid <- read.table(paste0(rt,".rsid"),col.names=c("snp","rsid"))
        suppressMessages(require(dplyr))
        suppressMessages(require(MendelianRandomization))
        if (file.exists(paste0(INF,"/mr/gsmr/trait/",p,"-",efo,"-top.txt")))
        {
          top <- read.table(paste0(INF,"/mr/gsmr/trait/",p,"-",efo,"-top.txt"),col.names=c("SNP","A1","A2","freq","b","se","p","N")) %>% pull(SNP)
          ld <- ieugwasr::ld_matrix(c(r,top),pop="EUR",with_alleles=FALSE)
          write.table(ld,file=paste0(rt,"-",r,".top"),quote=FALSE,sep="\t")
        }
        protein <- filter(gap.datasets::inf1,prot==p) %>% mutate(target.short=gsub("VEGF_A","VEGF-A",target.short)) %>% pull(target.short)
        source(file.path(INF,"rsid","gsmr_plot.r"))
        eff_dat <- paste0(rt,".eff_plot.gz")
        gsmr_data <- read_gsmr_data(eff_dat)
        gsmr_summary(gsmr_data)
        pdf(paste0(INF,"/mr/gsmr/out/",p,"-",trait,".eff_plot.pdf"))
        par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
        plot_gsmr_effect(gsmr_data, p, efo, colors()[75])
        dev.off()
        attach(gsmr_data)
        gsmr_snp_effect <- as.data.frame(snp_effect[,-(2:3)]) %>%
                           setNames(c("snp","effect_allele","other_allele","eaf","bx","bxse","by","byse")) %>%
                           mutate(bx=as.numeric(bx),bxse=as.numeric(bxse),by=as.numeric(by),byse=as.numeric(byse))
        j <- left_join(gsmr_snp_effect,snpid_rsid) %>% pull(rsid)
        rsid <- c(r,j)
        ldr <- ieugwasr::ld_matrix(rsid,pop="EUR",with_alleles=FALSE)
        write.table(ldr,file=paste0(rt,"-",r,".r"),quote=FALSE,sep="\t")
        detach(gsmr_data)
        attach(gsmr_snp_effect)
        png(paste0(INF,"/mr/gsmr/out/",p,"-",trait,".eff_plot.png"),res=300,height=7,width=7,units="in")
        mr_dat <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, exposure = protein, outcome = trait, snps = "snp",
                  effect_allele = effect_allele, other_allele = other_allele, eaf = eaf)
        mr_plot(mr_dat, interactive = FALSE, labels = FALSE)
        dev.off()
        detach(gsmr_snp_effect)
      END
    done
  done < <(sed '1d' ${EFO_UPDATE} | cut -f1,2,3,4)
}

function mr()
# now only a subset from OpenGWAS
{
  export prot=$(awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')
  export rt=${INF}/mr/gsmr
# gunzip -c ${rt}/prot/${prot}.gz | sed '1d' | cut -d' ' -f1 | sort -k1,1 | join - ${INF}/work/INTERVAL.rsid > ${rt}/prot/${prot}-rsid.txt
  Rscript -e '
    options(width=200)
    HPC_WORK <- Sys.getenv("HPC_WORK")
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rt <- Sys.getenv("rt")
#   ids <- scan(file.path(INF,"TNFB","efo"),"")
    pkgs <- c("TwoSampleMR","dplyr","ggplot2","openxlsx","stringr")
    invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))
    finngen <- c("finngen_R7_M13_ANKYLOSPON","finngen_R7_D3_SARCOIDOSIS","finngen_R7_M13_SJOGREN","finngen_R7_AB1_VIRAL_MENINGITIS","finngen_R7_AB1_HIV")
    gwas_catalog <- c("GCST90014325","GCST90134602","GCST90010715","GCST90061440","GCST90019016","GCST90014023")
    ids <- read.delim(file.path(INF,"OpenGWAS","efo-update.txt")) %>%
           pull(id) %>%
           setdiff(c(finngen,gwas_catalog))
    rsids <- read.table(file.path(rt,"prot",paste0(prot,"-rsid.txt")), col.names=c("snpid","SNP"))
    exposure_dat <- read_exposure_data(filename = file.path(rt,"prot",paste0(prot,".gz")), sep = " ",
      snp_col = "SNP",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "freq",
      beta_col = "b",
      se_col = "se",
      pval_col = "p",
      samplesize_col = "N"
    ) %>%
    mutate(region=gsub("chr|_|[a-z]","",SNP),snpid=gsub("(_[a-z])","\\U\\1",SNP,perl=TRUE)) %>%
    select(-SNP)
  # exposure_dat <- clump_data(left_join(exposure_dat,rsids))
    snpid_rsid <- extract_outcome_data(exposure_dat$region, ids) %>%
                  mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome)) %>%
                  select(snpid,SNP) %>%
                  distinct()
    exposure_dat <- clump_data(left_join(exposure_dat,snpid_rsid))
    outcome_dat <- extract_outcome_data(exposure_dat$region, ids) %>%
                   mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome),
                          effect_allele.outcome=toupper(effect_allele.outcome),
                          other_allele.outcome=toupper(other_allele.outcome)) %>%
                   select(snpid,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,
                          beta.outcome,se.outcome,pval.outcome,samplesize.outcome,id.outcome,outcome) %>%
                   group_by(id.outcome,SNP) %>%
                   slice(which.min(pval.outcome)) %>%
                   data.frame()
    dat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  # run_TwoSampleMR(harmonise, mr_plot="pQTLtools", prefix=prefix)
    mr_ivw_results <- mr(dat,method="mr_ivw") %>%
                      mutate(prot=prot,disease=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                      select(prot,id.outcome,b,se,pval,nsnp,disease)
    write.table(mr_ivw_results,file=file.path(rt,"mr",paste0(prot,".mr")),row.names=FALSE,quote=FALSE,sep="\t")
    mr_heterogeneity_results <- mr_heterogeneity(dat, method_list = "mr_egger_regression") %>%
                                mutate(prot=prot,disease=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                                select(prot,id.outcome,disease,method,Q,Q_df,Q_pval,-id.exposure,-exposure)
    write.table(mr_heterogeneity_results,file=file.path(rt,"mr",paste0(prot,".het")),row.names=FALSE,quote=FALSE,sep="\t")
    mr_single_results <- mr_singlesnp(dat) %>%
                      mutate(prot=prot,id.outcome=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                      select(prot,id.exposure,id.outcome,SNP,b,se,p)
    fp <- mr_forest_plot(mr_single_results)
    fp_outcome <- unlist(lapply(strsplit(names(fp),"[.]"),"[",2))
    for(i in 1:length(fp))
    {
      ggsave(fp[[i]],filename=file.path(rt,"mr",paste0(prot,"-",fp_outcome[i],".png")),device="png")
    }
  '
}

function mr2()
# This uses INTERVAL snpid/rsid pairings as well as OpenGWAS/GWAS Catalogu/FinnGen curations
{
  export prot=$(awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk '!/TNFB/' | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')
  export rt=${INF}/mr/gsmr
# gunzip -c ${rt}/prot/${prot}.gz | sed '1d' | cut -d' ' -f1 | sort -k1,1 | join - ${INF}/work/INTERVAL.rsid > ${rt}/prot/${prot}-rsid.txt
  Rscript -e '
    options(width=200)
    HPC_WORK <- Sys.getenv("HPC_WORK")
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rt <- Sys.getenv("rt")
    pkgs <- c("TwoSampleMR","dplyr","ggplot2","openxlsx","stringr")
    invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))
    gwas <- read.delim(file.path(INF,"OpenGWAS","efo-update.txt"))
    rsids <- read.table(file.path(rt,"prot",paste0(prot,"-rsid.txt")), col.names=c("SNP","rsid"))
    exposure_dat <- read_exposure_data(filename = file.path(rt,"prot",paste0(prot,".gz")), sep = " ",
                    snp_col = "SNP",
                    effect_allele_col = "A1",
                    other_allele_col = "A2",
                    eaf_col = "freq",
                    beta_col = "b",
                    se_col = "se",
                    pval_col = "p",
                    samplesize_col = "N") %>%
                    mutate(SNP=gsub("(_[a-z])","\\U\\1",SNP,perl=TRUE),id.exposure=prot,exposure=prot) %>%
                    left_join(rsids) %>%
                    mutate(SNP=rsid)
    user_dat <- function(gwas_id,gwas_trait)
    z <- read_outcome_data(filename=file.path(rt,"trait",paste0(prot,"-",gwas_id,"-rsid.txt")), sep=" ",
                           phenotype = gwas_id,
                           snp_col = "SNP",
                           beta_col = "b",
                           se_col = "se",
                           eaf_col = "freq",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "p",
                           samplesize_col = "N",
                           log_pval = FALSE
                          ) %>% mutate(id.outcome=gwas_id,outcome=gwas_trait)
    outcome_dat <- data.frame()
    for (gwas_id in pull(gwas,id)) outcome_dat <- bind_rows(outcome_dat,user_dat(gwas_id,filter(gwas,id==gwas_id)%>%pull(trait)))
    dat <- harmonise_data(clump_data(exposure_dat), outcome_dat, clump_r2 = 0.1, action=2)
  # run_TwoSampleMR(harmonise, mr_plot="pQTLtools", prefix=prefix)
    mr_ivw_results <- mr(dat,method="mr_ivw") %>%
                      mutate(prot=prot,disease=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                      select(prot,id.outcome,b,se,pval,nsnp,disease)
    write.table(mr_ivw_results,file=file.path(rt,"mr",paste0(prot,".mr")),row.names=FALSE,quote=FALSE,sep="\t")
    mr_heterogeneity_results <- mr_heterogeneity(dat, method_list = "mr_egger_regression") %>%
                                mutate(prot=prot,disease=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                                select(prot,id.outcome,disease,method,Q,Q_df,Q_pval,-id.exposure,-exposure)
    write.table(mr_heterogeneity_results,file=file.path(rt,"mr",paste0(prot,".het")),row.names=FALSE,quote=FALSE,sep="\t")
    mr_single_results <- mr_singlesnp(dat) %>%
                         mutate(prot=prot,id.outcome=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                         select(prot,id.exposure,id.outcome,SNP,b,se,p)
    write.table(mr_single_results,file=file.path(rt,"mr",paste0(prot,".single")),row.names=FALSE,quote=FALSE,sep="\t")
    fp <- mr_forest_plot(mr_single_results)
    fp_outcome <- gsub("^.","",gsub(prot,"",names(fp)))
    for(i in 1:length(fp))
    {
      ggsave(fp[[i]],filename=file.path(rt,"mr",paste0(prot,"-",fp_outcome[i],".png")),device="png")
    }
  '
}

function gsmr_trait()
{
# export prot=$(awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')
  export prot=$(sort -k4,4n ${INF}/TNFB/cis.dat  | awk '!/TNFB/' | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"] {print $3}')
  export region=$(sort -k4,4n ${INF}/TNFB/cis.dat | awk '!/TNFB/' | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"] {print $4":"$5"-"$6}')
  export rt=${INF}/mr/gsmr
  Rscript -e '
    options(width=200)
    HPC_WORK <- Sys.getenv("HPC_WORK")
    INF <- Sys.getenv("INF")
    prot <- Sys.getenv("prot")
    rt <- Sys.getenv("rt")
    region <- Sys.getenv("region")
#   ids <- scan(file.path(INF,"TNFB","efo"),"")
    suppressMessages(library(dplyr))
    ids <- read.delim(file.path(INF,"OpenGWAS","efo-update.txt")) %>% pull(opengwasid)
    pkgs <- c("TwoSampleMR","dplyr","ggplot2","openxlsx","stringr")
    invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))
    exposure_dat <- read.table(file.path(rt,"prot",paste0(prot,".gz")), header=TRUE) %>%
                    mutate(region=gsub("chr|_|[A-Z]","",SNP))
    outcome_dat <- extract_outcome_data(region, ids) %>%
                   mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome),
                          effect_allele.outcome=toupper(effect_allele.outcome),
                          other_allele.outcome=toupper(other_allele.outcome)
                         ) %>%
                   select(snpid,effect_allele.outcome,other_allele.outcome,eaf.outcome,
                          beta.outcome,se.outcome,pval.outcome,samplesize.outcome,id.outcome) %>%
                   setNames(c("SNP","A1","A2","freq","b","se","p","N","id.outcome")) %>%
                   group_by(id.outcome,SNP) %>%
                   slice(which.min(p)) %>%
                   data.frame()
    for (id in ids)
    {
      fp <- file.path(rt,"trait",paste0(prot,"-",id,".gz"))
      gz <- gzfile(fp,"w")
      write.table(filter(outcome_dat,id==id.outcome) %>% select(-id.outcome),gz,quote=FALSE,row.names=FALSE)
      cat(id,fp,file=file.path(rt,"trait",paste0("gsmr_",prot,"-",id)))
    }
  '
}
## EFO traits
# bfile
# mr2
# gsmr2
# gsmr_trait
gsmr3

function gsmr_hgi()
{
  for trait in A2 B1 B2 C2
  do
  export trait=${trait}
  export p=$(awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | \
             grep -w -f - ${INF}/work/INF1.merge.genes | \
             awk -vM=1e6 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $2}')
  if [ ! -f ${INF}/mr/gsmr/ref/gsmr_ref_${p} ]; then echo ${INF}/mr/gsmr/ref/${p} > ${INF}/mr/gsmr/ref/gsmr_ref_${p}; fi
  if [ ! -f ${INF}/mr/gsmr/prot/gsmr_${p} ]; then echo ${p} ${INF}/mr/gsmr/prot/${p}.gz > ${INF}/mr/gsmr/prot/gsmr_${p}; fi
  if [ ! -f ${INF}/mr/gsmr/trait/gsmr_${trait}_${p} ]; then
     echo ${trait} ${INF}/mr/gsmr/trait/${trait}-${p}.gz > ${INF}/mr/gsmr/trait/gsmr_${trait}_${p}
  fi
  gcta-1.9 --mbfile ${INF}/mr/gsmr/ref/gsmr_ref_${p} \
           --gsmr-file ${INF}/mr/gsmr/prot/gsmr_${p} ${INF}/mr/gsmr/trait/gsmr_${trait}_${p} --gsmr-direction 0 \
           --clump-r2 0.1 --gwas-thresh 5e-8 --diff-freq 0.4 --heidi-thresh 0.05 --gsmr-snp-min 10 --effect-plot \
           --out ${INF}/mr/gsmr/hgi//gsmr_${trait}_${p}
  if [ ! -f ${INF}/mr/gsmr/hgi/${trait}-${prot}.eff_plot.gz ]; then continue; fi
  R --no-save -q <<\ \ END
    INF <- Sys.getenv("INF")
    trait <- Sys.getenv("trait")
    p <- Sys.getenv("p")
    source(file.path(INF,"rsid","gsmr_plot.r"))
    gsmr_data <- read_gsmr_data(paste0(INF,"/mr/gsmr/hgi/gsmr_",trait,"_",p,".eff_plot.gz"))
    gsmr_summary(gsmr_data)
    pdf(paste0(INF,"/mr/gsmr/hgi/gsmr_",trait,"_",p,".eff_plot.pdf"))
    par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
    plot_gsmr_effect(gsmr_data, p, trait, colors()[75])
    dev.off()
  END
  done
}

# gsmr_hgi

### --- the following is obsolete

export p=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]{print $1}' ${INF}/work/inf1.tmp)

function MR_dat()
{
cut -f3 work/INF1.METAL | sed '1d' | sort | uniq | grep -w -f - work/INF1.merge.genes | awk -vjob=${job} 'NR==job' | \
parallel -j1 -C' ' '
  echo --- {2} ---
  gunzip -c METAL/{2}-1.tbl.gz | \
  cut -f1-6,10-12,18 | \
  awk -vchr={3} -vstart={4} -vend={5} -vM=1e6 -vlogp=-5.45131 -vsuffix=${suffix} "
        (suffix==\"cis\" && \$1==chr && \$2>=start-M && \$2 <= end+M && \$9<=logp) || (suffix==\"pan\" && \$9<=logp)
      " > work/mr/{2}-${suffix}.mri
  cut -f3 work/mr/{2}-${suffix}.mri > work/mr/{2}-${suffix}.mrs
  plink --bfile INTERVAL/cardio/INTERVAL --extract work/mr/{2}-${suffix}.mrs \
        --geno 0.1 --mind 0.1 --maf 0.005 --indep-pairwise 1000kb 1 0.01 --out work/mr/{2}-${suffix}
  (
    echo -e "rsid\tChromosome\tPosition\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tlogP\tN"
    grep -w -f work/mr/{2}-${suffix}.prune.in work/mr/{2}-${suffix}.mri | \
    awk "{\$3=\"chr\"\$1\":\"\$2;print}" | \
    sort -k3,3 | \
    join -23 -12 work/snp_pos - | \
    cut -d" " -f1 --complement | \
    tr " " "\t"
  ) | gzip -f > work/mr/{2}-${suffix}.mrx
'
}

function MR_dat_run()
{
  if [ ! -f work/INF1.merge.genes ]; then
     cut -f3,8,9,10 doc/olink.inf.panel.annot.tsv | grep -v BDNF | sed 's/"//g' | sort -k1,1 | join -12 work/inf1.tmp - > work/INF1.merge.genes
  fi
  if [ ! -d work/mr ]; then mkdir work/mr; fi
  for type in cis pan; do export suffix=${type}; MR_dat; done
}

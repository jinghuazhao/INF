suppressMessages(library(TwoSampleMR))
suppressMessages(library(dplyr))

# Single-SNP results

ao <- available_outcomes()
exposure_dat <- read_exposure_data(
 filename = 'TNFB.csv',
 sep = ',',
 snp_col = 'SNP',
 beta_col = 'beta',
 se_col = 'se',
 effect_allele_col = 'effect_allele',
 phenotype_col = 'Phenotype',
 units_col = 'units',
 other_allele_col = 'other_allele',
 eaf_col = 'eaf',
 samplesize_col = 'samplesize',
 ncase_col = 'ncase',
 ncontrol_col = 'ncontrol',
 gene_col = 'gene',
 pval_col = 'pval'
) %>%
mutate(SNP="rs2229092")
ids <- scan("efo","")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, ids,
               proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3) %>%
               distinct()
options(width=200)
exposure_dat
outcome_dat
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
save(dat,file="TNFB.rda")
mr_results <- mr(dat) %>%
              mutate(des=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
              select(id.outcome,b,se,pval,des)
write.table(mr_results,file="TNFB.txt",row.names=FALSE,sep="\t")

# IVW results

rt <- "~/INF/mr/gsmr"
prot <- "TNFB"

exposure_dat <- read_exposure_data(filename = file.path(rt,"prot",paste0(prot,".gz")), sep = ' ',
  snp_col = 'SNP',
  effect_allele_col = 'A1',
  other_allele_col = 'A2',
  eaf_col = 'freq',
  beta_col = 'b',
  se_col = 'se',
  pval_col = 'p',
  samplesize_col = 'N'
) %>%
mutate(region=gsub("chr|_|[a-z]","",SNP),snpid=gsub("(_[a-z])","\\U\\1",SNP,perl=TRUE)) %>%
select(-SNP)
outcome_dat <- extract_outcome_data(exposure_dat$region, ids,
               proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
save(outcome_dat,file="TNFB-outcome.rda")
outcome_dat <- outcome_dat %>%
               distinct() %>%
               mutate(snpid=gap::chr_pos_a1_a2(chr,pos,effect_allele.outcome,other_allele.outcome),SNP,
                      effect_allele.outcome=toupper(effect_allele.outcome),
                      other_allele.outcome=toupper(other_allele.outcome)) %>%
               select(snpid,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,
                      beta.outcome,se.outcome,pval.outcome,samplesize.outcome,id.outcome,outcome) %>%
               group_by(id.outcome,SNP) %>%
               slice(which.min(pval.outcome)) %>%
               data.frame()
exposure_dat <- left_join(exposure_dat,select(outcome_dat,snpid,SNP))

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
save(dat,file="TNFB-hamonise.rda")
mr_ivw_results <- mr(dat,method="mr_ivw") %>%
                  mutate(des=gsub("[ |]* id:[a-z]*-[a-z]-[A-Z]*[0-9]*|","",outcome)) %>%
                  select(id.outcome,b,se,pval,des)
write.table(mr_ivw_results,file="TNFB-ivw.txt",row.names=FALSE,quote=FALSE,sep="\t")

# not tested

f <- "TNFB.csv"
ivs <- format_data(read.csv(f))
library(pQTLtools)
pqtlMR(ivs, ids, mr_plot=FALSE)

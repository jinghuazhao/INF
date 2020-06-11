#!/usr/bin/bash

export pd=/rds/project/jmmh2/rds-jmmh2-public_databases/
export covid19=${pd}/covid19
export covid19sumstats=${pd}/covid19/COVID19_HGI_ANA5_20200429.txt
export getx7=${pd}/GTEx/GTEx_Analysis_v7_eQTL_all_associations

R --no-save -q <<END
  require(hyprcoloc)
  # Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)
  betas <- hyprcoloc::test.betas
  head(betas)
  ses <- hyprcoloc::test.ses
  head(ses)

  # Trait names and SNP IDs
  traits <- paste0("T", 1:10)
  rsid <- rownames(betas)

  # Colocalisation analyses
  results <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
END

zcat Whole_Blood.allpairs.txt.gz | cut -f1 | awk 'NR>1{print substr($1,1,15)}' | sort | uniq > work/wb
R --no-save -q <<END
   wb <- scan("work/wb", what="")
   g <- within(subset(grex::grex(wb),!is.na(uniprot_id)),{uniprot=trimws(uniprot_id)})
   m <- merge(gap::inf1,g,by="uniprot",all.x=TRUE)
END

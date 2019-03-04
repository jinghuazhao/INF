{ 
  OFS="\t"
  if (NR==1) print "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
  if ($18!="NA" && $18>=0.4) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $14
}

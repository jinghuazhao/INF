{ 
  OFS="\t"
  if (NR==1) print "SNPID", "CHR", "POS", "STRAND", "N", "EFFECT_ALLELE", "REFERENCE_ALLELE", "CODE_ALL_FQ", "BETA", "SE", "PVAL", "RSQ", "RSQ_IMP", "IMP"
<<<<<<< HEAD
  if ($18>=0.3) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $14
=======
  if ($18!="NA" && $18>=0.4) print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $14
>>>>>>> 34ec7abfb019ddd259c9ab69e517adc63c309bbd
}

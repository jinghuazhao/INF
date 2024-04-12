# CD40/IL6 forest plots

Lookup was done with PhenoScanner V2

```bash
module load ceuadmin/phenoscanner
phenoscanner --snp=rs2228145 -c All -x EUR -r 0.8
```

and crosscheck with GWAS Catalog, e.g.,

<https://www.ebi.ac.uk/gwas/publications/26482879>, under Variant and risk allele, we have rs2228145-C.

Two-sample MR as with forest plots are implmented in `10-4-2024.sh`.

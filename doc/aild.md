## Use of approximately independent LD blocks

The implementation involves several steps:

1. Set up 1703 autosomal regions as defined in [EURLD.bed](../tryggve/EURLD.bed).
2. Extract variants outside the [regions in high LD](../tryggve/high-LD-regions-hg19.txt) to [1672 regions](../tryggve/EURLD-no-high-LD-regions-hg19.bed) by [EURLD.sh](../tryggve/EURLD.sh).
3. Overlap regions and GWAS sumstats:
   * Tag GWAS sumstats with regions through [aild-rma.sb](../cardio/aild-rma.sb).
   * Pair protein-region which contains genomewide significant signals by [aild-list.sb](../cardio/aild-list.sb).
   * Independently, list variants by region in the reference panel by [aild-snplist.sb](../cardio/aild-snplist.sb).
4. clump via [aild-clump.sb](../cardio/aild-clump.sb).
5. cojo via [aild-cojo.sb](../cardio/aild-cojo.sb).
6. Downstream analyses with PhenoScanner (preferably v2) as in [ps.sh](../cardio/ps.sh) and forest plots on TRYGGVE (with study-specific sumstats).

The regions are predefined. As shown in [EURLD.tsv](../tryggve/EURLD.tsv) by [EURLD.R](../tryggve/EURLD.R), the LD patterns across the genome are more variable than the norm in a typical genomewide association analysis therefore slide windows such as 250kb (36), 500kb (300), or even 10Mb (1071), seeing that the sentinel variant may not necessarily lie right in the middle of a window. The number of signals in our case were close to GCTA but overestimated (53 by PLINK) as in [clump-cojo.md](../cardio/clump-cojo.md). For instance, it is often seen from the PLINK --clump-range outputs that sliding windows can give results in two neighbouring LD blocks.

Note that pairing regions of interest would reduce the burden of genomewide analysis, and also that region-specific reference will not affect results from steps 4 and 5 regarding use of variants from GWAS sumstats.

Steps 4 and 5 both use `INF1.aild`, which contains all the protein-region pairs. The results are classified as in [analysis.sh](../tryggve/analysis.sh). In particular, for step 5 this is done with [aild.sh](../cardio/aild.sh).

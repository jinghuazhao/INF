## Use of approximately independent LD blocks

The regions are predefined. As shown in [EURLD.tsv](tryggve/EURLD.tsv) by [EURLD.R](tryggve/EURLD.R), the LD patterns across the genome are more variable than the norm in a typical genomewide association analysis and it not reasonable to set slide windows such as 250kb (36), 500kb (300), or even 10Mb (1071), considering the sentinel variant may not necessarily lie right in the middle of a window. The number of signals would be underestimated (by GCGA) or overestimated (by PLINK) as in [cardio.md](cardio/cardio.md).

Steps to use:

1. Set up 1703 autosomal regions as defined in [EURLD.bed](tryggve/EURLD.bed).
2. Extract variants outside the [12 regions in high LD](tryggve/high-LD-regions-hg19.txt) to [1672 regions](tryggve/EURLD-no-high-LD-regions-hg19.bed) by [EURLD.sh](tryggve/EURLD.sh).
3. Overlap regions and GWAS sumstats:
   * Tag GWAS sumstats with regions through [aild-rma.sb](cardio/aild-rma.sb).
   * List variants by region in the reference panel by [aild-snplist.sb](cardio/aild-snplist.sb); to facilitate variant annotation only SNPs are kept.
   * Pair protein-region which contains genomewide significance by [aild-list.sb](cardio/aild-list.sb).
4. clump via [aild-clump.sb](cardio/aild-clump.sb).
5. cojo via [aild-cojo.sb](cardio/aild-cojo.sb).

Results from 4 and 5 are classified as in [analysis.sh](tryggve/analysis.sh).

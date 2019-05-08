#!/bin/bash

awk -vOFS="\t" '{
  if(NR==1) print "#chrom","Start","End", "r"
  print "chr" $1, $2, $3, $4
}' tryggve/high-LD-regions-hg19.txt > tryggve/high-LD-regions-hg19.bed

module load gcc/5.2.0

(
  head -1 tryggve/EURLD.bed
  bedtools intersect -v -a tryggve/EURLD.bed -b tryggve/high-LD-regions-hg19.bed
) > tryggve/EURLD-no-high-LD-regions-hg19.bed

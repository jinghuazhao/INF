#!/bin/bash

module load gcc/5.2.0

cd tryggve

awk -vOFS="\t" '{
  if(NR==1) print "#chrom","Start","End", "r"
  print "chr" $1, $2, $3, $4
}' high-LD-regions-hg19.txt > high-LD-regions-hg19.bed

(
  head -1 EURLD.bed
  bedtools intersect -v -a EURLD.bed -b high-LD-regions-hg19.bed
) > EURLD-no-high-LD-regions-hg19.bed

cd -

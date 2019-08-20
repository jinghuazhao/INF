{
   d3=$13; gsub(/?/,"",d3)
   if (length(d3) >= 3 && $18 >= 3500)
   if ($12 > -9.30103) print;
   else {
      if ($14 < 30) print;
      else {
        d3n=d3; d3p=d3;
        gsub(/+|-|p/,"",d3n); gsub(/+|-|n/,"",d3p);
        if (length(d3n) >= 3 || length(d3p) >= 3) print;
      }
   }
}
# R
# > log10(5e-10)
# [1] -9.30103
# head -1 METAL/4E.BP1-1.tbl | sed 's|\t|\n|g' | awk '{print "#" NR,$1}'
#1 Chromosome
#2 Position
#3 MarkerName
#4 Allele1
#5 Allele2
#6 Freq1
#7 FreqSE
#8 MinFreq
#9 MaxFreq
#10 Effect
#11 StdErr
#12 log(P)
#13 Direction
#14 HetISq
#15 HetChiSq
#16 HetDf
#17 logHetP
#18 N

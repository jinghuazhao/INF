{
   d3=$13; gsub(/?/,"",d3); l=length(d3); s=sprintf("%*s", l, "")
   d3r=d3
   gsub(/n/,"-",d3r); n=gensub(/ /, "-", "g", s)
   gsub(/p/,"+",d3r); p=gensub(/ /, "+", "g", s)
   if (l >= 3 && $18 >= 3500)
      if ($12 > -9.30103) print;
      else {
         if ($14 < 30) print;
         else if (d3 == "nnn" || d3 == "ppp") print;
         else if (d3r == n || d3r == p ) print
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

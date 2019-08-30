# 29-8-2019 JHZ

if [ -f work/INF1.gsmr ]; then rm work/INF1.gsmr; fi
(
  cat work/gsmr_CAD*.gsmr | \
  head -1
  ls work/gsmr_CAD*gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > INF1_CAD.gsmr

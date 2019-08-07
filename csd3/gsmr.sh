# 7-8-2019 JHZ

if [ -f work/INF1.gsmr ]; then rm work/INF1.gsmr; fi
(
  cat work/*.gsmr | \
  head -1
  ls work/*.gsmr | \
  parallel -j1 -C' ' '
    if [ -f {} ]; then
       awk "NR>1" {}
    fi
  '
) | \
grep -v nan > INF1.gsmr

# 17-7-2019 JHZ

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

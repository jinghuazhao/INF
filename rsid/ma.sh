#!/usr/bin/bash

# snpid --> rsid
function primary()
{
  rm -f ${INF}/sentinels/*-rsid.ma
  ls ${INF}/sentinels/*.ma | \
  xargs -I {} sh -c "basename {} .ma" | \
  parallel --env INF -C' ' '
  (
    head -1 ${INF}/sentinels/{1}.ma
    sed "1d" ${INF}/sentinels/{1}.ma | \
    sort -k1,1 | \
    join ${INF}/work/INTERVAL.rsid - | \
    awk "{\$1=\"\";print}" | \
    awk "{\$1=\$1};1"
  ) > ${INF}/sentinels/{1}-rsid.ma
  '
}

function all()
{
  rm -f ${INF}/work/*-rsid.ma
  ls ${INF}/work/*.ma | \
  xargs -I {} sh -c "basename {} .ma" | \
  parallel --env INF -C' ' '
  (
    head -1 ${INF}/work/{1}.ma
    sed "1d" ${INF}/work/{1}.ma | \
    sort -k1,1 | \
    join ${INF}/work/INTERVAL.rsid - | \
    awk "{\$1=\"\";print}" | \
    awk "{\$1=\$1};1"
  ) > ${INF}/work/{1}-rsid.ma
  '
}

cd work
primary
cd -

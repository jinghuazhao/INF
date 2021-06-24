#!/usr/bin/bash

# snpid --> rsid
function primary()
{
  ls ${INF}/METAL/*-1.tbl.gz | \
  xargs -I {} sh -c "basename {} -1.tbl.gz" | \
  parallel --env INF -C' ' '
    export p={}
    echo {}
    (
      echo SNP A1 A2 freq b se p N
      zcat ${INF}/METAL/{}-1.tbl.gz | \
      awk "10^\$12 <=5e-10 {print \$3,toupper(\$4),toupper(\$5),\$6,\$10,\$11,10^\$12,\$18}" | \
      sort -k1,1 | \
      join ${INF}/work/INTERVAL.rsid - | \
      awk "{\$1=\"\";print}" | \
      awk "{\$1=\$1};1"
    ) > ${INF}/sentinels/{}-rsid.ma
  '
}

function all()
{
  rm -f ${INF}/work/*-rsid.ma
  ls ${INF}/work/*.ma | \
  xargs -I {} sh -c "basename {} .ma" | \
  parallel --env INF -C' ' '
  (
    head -1 ${INF}/work/{}.ma
    sed "1d" ${INF}/work/{}.ma | \
    sort -k1,1 | \
    join ${INF}/work/INTERVAL.rsid - | \
    awk "{\$1=\"\";print}" | \
    awk "{\$1=\$1};1"
  ) > ${INF}/work/{}-rsid.ma
  '
}

primary

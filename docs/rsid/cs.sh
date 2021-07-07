#!/usr/bin/bash

if [ ! -d ${INF}/cs ]; then mkdir ${INF}/cs; fi

function tbl()
# all variants in the region
{
  cut -f5,6,8,9 --output-delimiter=' ' ${INF}/work/INF1.merge | sed '1d' | \
  parallel -C' ' --env INF '
    echo {1}-{2}-{3}-{4}
    export p={1}
    export r={2}
    export chr={3}
    export pos={4}
    zcat ${INF}/METAL/{1}-1.tbl.gz | \
    cut -f1-3,6,10-12,18 | \
    awk -vchr=${chr} -vpos=${pos} "NR==1 || (\$1==chr && \$2 >= pos - 1e6 && \$2 < pos + 1e6)" | \
    gzip -f > ${INF}/cs/{1}-{2}.tbl.gz
  R --no-save <<\ \ END
    options(width=200)
    p <- Sys.getenv("p")
    r <- Sys.getenv("r")
    pr <- paste0(p,"-",r)
    tbl <- read.delim(paste0(pr,".z"),sep=" ")
    z <- suppressMessages(library(gap))
    z <- suppressMessages(cs(tbl))
    write(z[["MarkerName"]],file=paste0(pr,".metal"),nrow(z))
  END
 '
}

cd ${INF}/cs
(
  awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
  parallel --env INF -C' ' '
      ${INF}/rsid/cs.inc {3}
      awk -vprot={1} -vrsid={2} -vOFS="\t" "{print prot,rsid,\$0}" ${INF}/cs/{1}-{2}.cs
  '
) > ${INF}/work/INF1.merge-rsid.cs

R --no-save -q <<END
  library(dplyr)
  INF <- Sys.getenv("INF")
  cs <- read.table(file.path(INF,"work","INF1.merge-rsid.cs"), col.names=c("prot","pQTL","set"), sep='\t')
  for (i in 1:nrow(cs))
  {
    cs$n[i] <- length(unlist(strsplit(cs$set[i]," ")))
    cs$isin[i] <- grepl(cs$pQTL[i],cs$set[i])
  }
  options(width=200)
  table(cs$n)
  table(cs$isin)
END
cd -

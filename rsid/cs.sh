#!/usr/bin/bash

if [ ! -d ${INF}/cs ]; then mkdir ${INF}/cs; fi

function tbl()
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
    p <- Sys.getenv("p")
    r <- Sys.getenv("r")
    pr <- paste0(p,"-",r)
    INF <- Sys.getenv("INF")
    tbl <- subset(read.delim(file.path(INF,"cs",paste0(pr,".tbl.gz"))),!is.na(log.P.))
    flag <- with(tbl,Freq1>0.5)
    tbl[flag,"Freq1"] <- 1-tbl[flag,"Freq1"]
    require(corrcoverage)
    ppa <- with(tbl,z0_pp(z=Effect/StdErr,f=Freq1,type="quant",N=N,s=0))
    idx <- with(credset(ppa, thr=0.95),credset)
    write.table(tbl[["MarkerName"]][idx],file=file.path(INF,"cs",paste0(pr,".cs")),col.names=FALSE,row.names=FALSE,quote=FALSE)
  END
 '
}

awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
parallel --env INF -C' ' '
  echo {1}-{2} {3}
  ${INF}/rsid/cs.inc {3}
'

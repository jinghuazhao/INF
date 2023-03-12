#!/usr/bin/bash

if [ ! -d ${INF}/cs ]; then mkdir ${INF}/cs; fi

function cs()
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
    awk -vchr=${chr} -vpos=${pos} "NR==1 || (\$1==chr && \$2 >= pos - 1e6 && \$2 < pos + 1e6) && !(chr==19 && pos >= 53296855 && pos <= 54500000)" | \
    gzip -f > ${INF}/cs/unprune/{1}-{2}.tbl.gz
  R --no-save <<\ \ END
    options(width=200)
    p <- Sys.getenv("p")
    r <- Sys.getenv("r")
    pr <- paste0(p,"-",r)
    tbl <- read.delim(paste0("~/INF/cs/unprune/",pr,".tbl.gz"))
    z <- suppressMessages(gap::cs(tbl))
    write(z[["MarkerName"]],file=paste0("~/INF/cs/unprune/",pr,".snpid"),nrow(z))
    write(z[["ppa"]],file=paste0("~/INF/cs/unprune/",pr,".ppa"),nrow(z))
  END
  cat ~/INF/cs/unprune/${p}-${r}.snpid | \
  tr " " "\n" | \
  sort -k1,1 | \
  join - ${INF}/work/INTERVAL.rsid | \
  cut -d" " -f2 | \
  tr "\n" " " > ~/INF/cs/unprune/${p}-${r}.cs
  '
  (
    awk 'NR > 1 {print $3,$2,$1}' ${INF}/work/INF1.METAL | \
    parallel --env INF -C' ' 'awk -vprot={1} -vrsid={2} -vOFS="\t" "{print prot,rsid,\$0}" ${INF}/cs/unprune/{1}-{3}.cs'
  ) > ${INF}/cs/unprune/INF1.merge-rsid.cs
  (
    awk 'NR > 1 {print $3,$2,$1}' ${INF}/work/INF1.METAL | \
    parallel --env INF -C' ' '
      cat ${INF}/cs/unprune/{1}-{3}.ppa | \
      Rscript -e "options(width=2000);ppa <- scan(\"stdin\");cat(format(ppa,digits=3,scientific=TRUE))" | \
      awk -vprot={1} -vrsid={2} -vOFS="\t" "{print prot,rsid,\$0}"
    '
  ) > ${INF}/cs/unprune/INF1.merge-rsid.ppa
}

function prune()
{
(
  awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
  parallel --env INF -C' ' '
      ${INF}/rsid/cs.inc {3}
      awk -vprot={1} -vrsid={2} -vOFS="\t" "{print prot,rsid,\$0}" ${INF}/cs/{1}-{2}.cs
  '
) > ${INF}/work/INF1.merge-rsid.cs
(
  awk 'NR > 1 {print $5,$6,NR-1}' ${INF}/work/INF1.merge-rsid | \
  parallel --env INF -C' ' '
      ${INF}/rsid/cs.inc {3}
      awk -vprot={1} -vrsid={2} -vOFS="\t" "{print prot,rsid,\$0}" ${INF}/cs/{1}-{2}.ppa
  '
) > ${INF}/work/INF1.merge-rsid.ppa
}

R --no-save -q <<END
  library(dplyr)
  INF <- Sys.getenv("INF")
  sumcs <- function(cs)
  {
    for (i in 1:nrow(cs))
    {
      cs$n[i] <- length(unlist(strsplit(cs$set[i]," ")))
      cs$isin[i] <- grepl(cs$pQTL[i],cs$set[i])
    }
    n <- table(cs$n)
    isin <- table(cs$isin)
    print(n)
    print(isin)
    list(n=n,isin=isin)
  }
  options(width=200)
  cat("prune\n")
  prune <- read.table(file.path(INF,"cs","prune","INF1.merge-rsid.cs"), col.names=c("prot","pQTL","set"), sep='\t')
  summary_prune <- sumcs(prune)
  cat("unprune\n")
  unprune <- read.table(file.path(INF,"cs","unprune","INF1.merge-rsid.cs"), col.names=c("prot","pQTL","set"), sep='\t')
  summary_unprune <- sumcs(unprune)
END

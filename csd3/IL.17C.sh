#!/usr/bin/bash

function setup()
# IL.17C
{
tabix sumstats/INTERVAL/INTERVAL.IL.17C.gz 16 | sort -k11,11gr | awk '$11!="NA"' | tail
(
cat <<END
chr16:88684495_G_T rs17700884
chr16:88695338_C_G rs4782388
chr16:88692206_A_ACGC 16:88692206_A_ACGC
chr16:88692204_T_TGAC rs71391399
chr16:88692203_T_TTGACGA rs559941365
END
) | \
sed 's/chr//;s/:/ /;s/_/ /g' | \
awk '{print $1,$2,$2,$3"/"$4,$5}' | \
parallel -C' ' '
vep --id "{1} {2} {3} {4}" -o {5} --species homo_sapiens --assembly GRCh37 --cache --offline --force_overwrite --tab --nearest symbol --pick
'
gunzip -c ${INF}/sumstats/INTERVAL/INTERVAL.IL.17C.gz | cut -f1-3,9-10 | gzip -f > ${INF}/work/IL.17C.gz

Rscript -e '
  suppressMessages(library(dplyr))
  suppressMessages(library(gap))
  INF <- Sys.getenv("INF")
  genes <- data.frame(MarkerName=c("chr16:88692203_T_TTGACGA"),gene=c("IL17C"),color=c("red"))
  gz <- gzfile(file.path(INF,"work","IL.17C.gz"))
  IL.17C <- read.delim(gz,as.is=TRUE) %>%
            mutate(Z=BETA/SE,P=pvalue(Z),log10P=-log10p(Z)) %>%
            rename(Chromosome=CHR,Position=POS,MarkerName=SNPID) %>%
            select(Chromosome,Position,MarkerName,Z,P,log10P) %>%
            left_join(genes)
  save(IL.17C, genes, file=file.path(INF,"work","IL.17C.rda"))
'

R --no-save -q <<END
# 24/9/2020, not highlighted
  library(gap)
# library(Rmpfr)
  INF <- Sys.getenv("INF")
  gz <- gzfile(file.path(INF,"METAL","IL.17C-1.tbl.gz"))
  IL.17C <- within(read.delim(gz,as.is=TRUE), {
   Z <- Effect/StdErr;
#  P <- as.numeric(2*pnorm(mpfr(abs(Z),100),lower.tail=FALSE))
   P <- 2*pnorm(abs(Z),lower.tail=FALSE)
  })
  subset(IL.17C, P==0)
  png("IL.17C.png", res=300, units="in", width=9, height=6)
  par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
  mhtplot.trunc(IL.17C, chr="Chromosome", bp="Position", p="P", snp="MarkerName", z = "Z",
                suggestiveline=FALSE, genomewideline=-log10(5e-10), logp = TRUE,
                cex.mtext=0.6, cex.text=0.7,
                mtext.line=4, y.brk1=30, y.brk2=50, cex.axis=1.2, cex=0.5, y.ax.space=10,
                col = c("blue4", "skyblue")
  )
  dev.off()
END
}

Rscript -e '
  suppressMessages(library(dplyr))
  suppressMessages(library(gap))
  INF <- Sys.getenv("INF")
  load(file.path(INF,"work","IL.17C.rda"))
  subset(IL.17C,!is.na(gene))
  png("IL.17C-mhtplot.trunc.png", res=300, units="in", width=9, height=6)
  par(oma=c(0,0,0,0), mar=c(5,6.5,1,1))
  mhtplot.trunc(subset(IL.17C,!is.na(Z),select=-color), chr="Chromosome", bp="Position", z="Z", snp="MarkerName",
                suggestiveline=-log10(1e-7), genomewideline=-log10(5e-10),
                cex.mtext=1.2, cex.text=0.7,
                annotatelog10P=-log10(6.35e-9), annotateTop = FALSE, highlight=with(genes,gene),
                mtext.line=3, y.brk1=2, y.brk2=3, trunc.yaxis=FALSE, delta=0.1, cex.axis=1.2,
                cex=0.5, font=2, font.axis=1, y.ax.space=1,
                col = c("blue4", "skyblue")
  )
  dev.off()
  mhtdata <- filter(IL.17C,!is.na(Z)) %>%
             select(-MarkerName,-Z,-log10P) %>%
             mutate(P=as.numeric(P))
  cis <- with(mhtdata,Chromosome==16 & Position>=88704999-1e6 & Position<=88706881+1e6)
  mhtdata[cis,"color"] <- seq(length(colors()))[colors()=="red"]
  subset(mhtdata,!is.na(gene))
  png("IL.17C-mhtplot.png", res=300, units="in", width=9, height=6)
  opar <- par()
  par(cex=0.5)
  ops <- mht.control(colors=rep(c("blue4","skyblue"),11),srt=0,yline=2.5,xline=2)
  hops <- hmht.control(data=subset(mhtdata,!is.na(gene)))
  mhtplot2(mhtdata,ops,hops,xlab="",ylab="",srt=0, cex.axis=1.2)
  axis(2,at=1:10)
  title("")
  par(opar)
  dev.off()
'

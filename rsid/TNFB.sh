#!/usr/bin/bash

module load gcc/6 texlive
Rscript -e "knitr::knit('TNFB.Rmd')"
pandoc TNFB.md --citeproc --mathjax -o TNFB.docx
pandoc TNFB.md --citeproc --mathjax --self-contained -o TNFB.html
scp TNFB.docx TNFB.html ivw-*png jhz22@shell.srcf.net:/home/jhz22/public_html/INF/latest/

function pdf()
{
  pandoc TNFB.md -f markdown --variable fontsize:12pt --variable=geometry:"paperwidth=18in, paperheight=12in, margin=24pt" -o TNFB.pdf
  scp TNFB.docx TNFB.html TNFB.pdf *png jhz22@shell.srcf.net:/home/jhz22/public_html/INF/latest/
}


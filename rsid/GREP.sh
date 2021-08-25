#!/usr/bin/bash

source ~/COVID-19/py37/bin/activate

cut -d, -f10,14 ~/INF/work/INF1.merge.cis.vs.trans | \
grep -v trans | \
cut -d, -f1 > inf1.genes

python grep.py --genelist inf1.genes --out inf1  --test ATC --output-drug-name

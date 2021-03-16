#!/usr/bin/bash

cd ${INF}/GREAT
cat <(head -1 IL12B.tsv) \
    <(sed -e 's/Ontolog.*/IL12B/' IL12B.tsv) \
    <(sed -e 's/Ontology.*/KITLG/' KITLG.tsv) \
    <(sed -e 's/Ontology.*/TNFSF10/' TNFSF10.tsv) > IL12B-KITLG-TNFSF10.tsv
cd -

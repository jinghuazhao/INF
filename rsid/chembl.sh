#!/usr/bin/bash

source ~/COVID-19/py37/bin/activate

export chembl=INF1.merge.chembl
export f=INF1.merge.compounds

python <<END
from chembl_webresource_client.new_client import new_client
import csv
import json
import os

# Get many compounds by their ChEMBL IDs:
chembl = os.environ['chembl']
chembl_ids = []
i = open(chembl, 'r')
with i:
    reader = csv.DictReader(i,delimiter='\t')
    for row in reader:
        chembl_id = row['To']
        chembl_ids.append(chembl_id)

print(chembl_ids)
f = os.environ['f']
molecule = new_client.molecule
records = molecule.get(chembl_ids)
with open(f, 'w') as outfile:
    json.dump(records, outfile)

# Search target by gene name:
target = new_client.target
gene_name = 'IL12B'
res = target.search(gene_name)
END

R --no-save <<END
require(jsonlite)
f <- Sys.getenv('f')
r <- read_json(f, simplifyVector = TRUE)
# write.table(r,file=paste0(f,"-txt"),quote=FALSE,row.names=FALSE,sep="\t")
END

# pip install chembl_webresource_client

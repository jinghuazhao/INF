#!/usr/bin/bash

source ~/COVID-19/py37/bin/activate
# pip install chembl_webresource_client

export f=INF1.merge.compounds

python <<END
from chembl_webresource_client.new_client import new_client
import json
import os

# Search target by gene name:
target = new_client.target
gene_name = 'IL12B'
res = target.search(gene_name)

# Get many compounds by their ChEMBL IDs:
molecule = new_client.molecule
records = molecule.get(['CHEMBL6498', 'CHEMBL6499', 'CHEMBL6505'])
f = os.environ['f']
with open(f, 'w') as outfile:
    json.dump(records, outfile)

END

R --no-save <<END
require(jsonlite)
f <- Sys.getenv('f')
r <- read_json(f, simplifyVector = TRUE)
# write.table(r,file=paste0(f,"-txt"),quote=FALSE,row.names=FALSE,sep="\t")
END

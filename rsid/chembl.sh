#!/usr/bin/bash

source ~/COVID-19/py37/bin/activate

export chembl=INF1.merge.chembl
export f=INF1.merge.compounds

(
  echo gene
  cut -d, -f10 INF1.merge.cis.vs.trans | sed '1d' | sort | uniq
) > INF1.merge.glist

python <<END
from chembl_webresource_client.new_client import new_client
import csv
import json
import os

# Search target by gene name:
target = new_client.target
gene_ids = []
res = []
g = open('INF1.merge.glist', 'r')
with g:
    reader = csv.DictReader(g)
    for row in reader:
        gene_id = row['gene']
        res.append(target.search(gene_id))
        gene_ids.append(gene_id)

print(gene_ids)
g.close()

from django.db.models.loading import get_model
def dump(qs, outfile_path):
	"""
	Takes in a Django queryset and spits out a CSV file.
	
	Usage::
	
		>> from utils import dump2csv
		>> from dummy_app.models import *
		>> qs = DummyModel.objects.all()
		>> dump2csv.dump(qs, './data/dump.csv')
	
	Based on a snippet by zbyte64::
		
		http://www.djangosnippets.org/snippets/790/
	
	"""
        model = qs.model
	writer = csv.writer(open(outfile_path, 'w'))
	
	headers = []
	for field in model._meta.fields:
		headers.append(field.name)
	writer.writerow(headers)
	
	for obj in qs:
		row = []
		for field in headers:
			val = getattr(obj, field)
			if callable(val):
				val = val()
			if type(val) == unicode:
				val = val.encode("utf-8")
			row.append(val)
		writer.writerow(row)

# Get targets by their ChEMBL IDs:
chembl = os.environ['chembl']
chembl_ids = []
i = open(chembl, 'r')
with i:
    reader = csv.DictReader(i,delimiter='\t')
    for row in reader:
        chembl_id = row['To']
        chembl_ids.append(chembl_id)

print(chembl_ids)
i.close()

f = os.environ['f']
target = new_client.target
records = target.get(chembl_ids)

# this requires R/jsonlite
with open('INF1.merge.target', 'w') as outfile:
    json.dump(records, outfile)

outfile.close()

data_file = open('INF1.merge.target', 'w')
csv_writer = csv.writer(data_file)
count = 0
for r in records:
    if count==0:
       header=r.keys()
       csv_writer.writerow(header)
       count += 1
    csv_writer.writerow(r.values())

data_file.close()
END

R --no-save <<END
  require(jsonlite)
  f <- Sys.getenv('f')
  r <- read_json(f, simplifyVector = TRUE)
  # write.table(r,file=paste0(f,"-txt"),quote=FALSE,row.names=FALSE,sep="\t")
END

# pip install chembl_webresource_client

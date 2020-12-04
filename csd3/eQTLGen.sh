export eQTLGen=~/rds/public_databases/eQTLgen

function cis_pQTL_cis_eQTL()
{
  grep -f work/INF1.merge.cis work/INTERVAL.rsid > work/INF1.merge.cis-rsid
(
python3 <<END
#!/usr/bin/python3
import csv 

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'])
END
) | \
sort -k2,2 | \
join -12 - <(sort -k1,1 work/INF1.merge.cis-rsid) | \
sort -k4,4 | \
join -14 -22 - <(cut -d' ' -f2 work/INF1.merge.cis-rsid | zgrep -f - -w ${eQTLGen}/cis.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
awk '$4==$8' > work/eQTLGen.cis
}

function trans_pQTL_trans_eQTL()
{
  cut -f6 work/INF1.merge.trans | grep -f - work/INTERVAL.rsid > work/INF1.merge.trans-rsid
(
python3 <<END
#!/usr/bin/python3
import csv

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'])
END
) | \
sort -k2,2 | \
join -12 - <(sort -k1,1 work/INF1.merge.trans-rsid) | \
sort -k4,4 | \
join -14 -22 - <(cut -d' ' -f2 work/INF1.merge.trans-rsid | zgrep -f - -w ${eQTLGen}/trans.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
awk '$4==$8' > work/eQTLGen.trans
}

function cis_eQTL()
{
(
python3 <<END
#!/usr/bin/python3
import csv

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'], row['cis.trans'])
END
) | \
grep cis | \
sort -k3,3 | \
join -13 -25 - <(gunzip -c ${eQTLGen}/cis.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > work/eQTLGen.cis-all
}

function trans_eQTL()
{
(
python3 <<END
#!/usr/bin/python3
import csv

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'], row['cis.trans'])
END
) | \
grep trans | \
sort -k3,3 | \
join -13 -25 - <(gunzip -c ${eQTLGen}/trans.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > work/eQTLGen.trans-all
}

R --no-save <<END
  ys1 <- c(paste0("Yes-",1:22),paste0("No2-",1:437))
  ys2 <- c(paste0("Yes-",1:22),paste0("No1-",1:155))
  ys <- list(ys1,ys2)
  names(ys) <- c("eQTL","pQTL")
  VennDiagram::venn.diagram(x = ys, filename='eQTLGen.png', imagetype="png", output=TRUE,
                            height=12, width=12, units="cm", resolution=500,
                            fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
END

export eQTLGen=~/rds/public_databases/eQTLgen

function cistrans_python()
{
python3 <<END
#!/usr/bin/python3
import csv 

with open('work/INF1.merge.cis.vs.trans', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['uniprot'], row['SNP'], row['p.gene'])
END
}

for cistrans in cis trans
do
  export cistrans_rsid=work/INF1.${cistrans}-rsid
  grep -f work/INF1.${cistrans} work/INTERVAL.rsid > ${cistrans_rsid}
  cistrans_python | \
  sort -k2,2 | \
  join -12 - <(sort -k1,1 ${cistrans_rsid}) | \
  sort -k4,4 | \
  join -14 -22 - <(cut -d' ' -f2 ${cistrans_rsid} | zgrep -f - -w ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
  awk '$4==$8' > work/eQTLGen.${cistrans}
  cistrans_python | \
  sort -k3,3 | \
  join -13 -25 - <(gunzip -c ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > work/eQTLGen.${cistrans}-all
done

R --no-save <<END
  ys1 <- c(paste0("Yes-",1:22),paste0("No1-",1:62238))
  ys2 <- c(paste0("Yes-",1:22),paste0("No2-",1:37))
  ys <- list(ys1,ys2)
  names(ys) <- c("eQTL","pQTL")
  VennDiagram::venn.diagram(x = ys, filename='work/eQTLGen-cis.png', imagetype="png", output=TRUE,
                            height=12, width=12, units="cm", resolution=500,
                            fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
  ys1 <- c(paste0("Yes-",1:3),paste0("No1-",1:1129))
  ys2 <- c(paste0("Yes-",1:3),paste0("No2-",1:118))
  ys <- list(ys1,ys2)
  names(ys) <- c("eQTL","pQTL")
  VennDiagram::venn.diagram(x = ys, filename='work/eQTLGen-trans.png', imagetype="png", output=TRUE,
                            height=12, width=12, units="cm", resolution=500,
                            fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
END

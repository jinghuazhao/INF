#!/usr/bin/bash

export eQTLGen=~/rds/public_databases/eQTLGen

function cistrans_python()
# select named columns from .csv
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
# all from eQTLGen but with LD
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

function jma()
{
if [ ! -d ${INF}/eQTLGen ]; then mkdir ${INF}/eQTLGen; fi

cd ${INF}
module load ceuadmin/stata 
stata <<END
  insheet using "sentinels/INF1.jma-rsid.cis.vs.trans.", case clear delim(" ")
  l uniprot SNP pgene
  outsheet uniprot SNP pgene using "eQTLGen/uniprot-rsid-gene.", delim(" ") noquote noname replace
  outsheet SNP using "eQTLGen/INF1.cis-rsid" if cistrans=="cis", noname noquote replace
  outsheet SNP using "eQTLGen/INF1.trans-rsid" if cistrans=="trans", noname noquote replace
END
cd -

for cistrans in cis trans
do
  export cistrans_rsid=${INF}/eQTLGen/INF1.${cistrans}-rsid
# cistrans_stata
  sort -k2,2 ${INF}/eQTLGen/uniprot-rsid-gene | \
  join -12 - <(sort -k1,1 ${cistrans_rsid}) | \
  sort -k1,1 | \
  join -11 -22 - <(cut -d' ' -f2 ${cistrans_rsid} | zgrep -f - -w ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k2,2) | \
  awk '$3==$7' > ${INF}/eQTLGen/eQTLGen.${cistrans}
# all from eQTLGen but with LD
  sort -k3,3 ${INF}/eQTLGen/uniprot-rsid-gene | \
  join -13 -25 - <(gunzip -c ${eQTLGen}/${cistrans}.txt.gz | cut -f1,2,5,6,9 | sort -k5,5) > ${INF}/eQTLGen/eQTLGen.${cistrans}-all
done

R --no-save -q <<END
  ys1 <- c(paste0("Yes-",1:37),paste0("No1-",1:81931))
  ys2 <- c(paste0("Yes-",1:37),paste0("No2-",1:62))
  ys <- list(ys1,ys2)
  names(ys) <- c("eQTL","pQTL")
  INF <- Sys.getenv("INF")
  VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"eQTLGen","eQTLGen-cis.png"), imagetype="png", output=TRUE,
                            height=12, width=12, units="cm", resolution=500,
                            fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
  ys1 <- c(paste0("Yes-",1:14),paste0("No1-",1:1462))
  ys2 <- c(paste0("Yes-",1:14),paste0("No2-",1:114))
  ys <- list(ys1,ys2)
  names(ys) <- c("eQTL","pQTL")
  VennDiagram::venn.diagram(x = ys, filename=file.path(INF,"eQTLGen","eQTLGen-trans.png"), imagetype="png", output=TRUE,
                            height=12, width=12, units="cm", resolution=500,
                            fill=c("yellow","purple"), cat.pos=c(-30,30), rotation.degree = 0)
END
}

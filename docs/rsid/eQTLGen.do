local inf : env INF
insheet using "`inf'/sentinels/INF1.jma-rsid.cis.vs.trans.", case clear delim(" ")
l uniprot SNP pgene
outsheet uniprot SNP pgene using "`inf'/eQTLGen/uniprot-rsid-gene.", delim(" ") noquote noname replace
outsheet SNP using "`inf'/eQTLGen/INF1.cis-rsid" if cistrans=="cis", noname noquote replace
outsheet SNP using "`inf'/eQTLGen/INF1.trans-rsid" if cistrans=="trans", noname noquote replace

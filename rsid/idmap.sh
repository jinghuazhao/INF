#!/usr/bin/bash

export gwas_id=${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_9.fam
export inf_id=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export scallop_inf_id=${HOME}/COVID-19/SCALLOP-Seq

wc -l ${gwas_id}
wc -l <(sed '1,2d' ${inf_id})
join <(cut -f1 ${gwas_id} | sort -k1,1n) <(sed '1,2d' ${inf_id} | cut -d' ' -f1 | sort -k1,1n) | wc -l

module load ceuadmin/stata

stata <<END
local scallop_inf_id : env scallop_inf_id
di "\`scallop_inf_id\'"
insheet using \`scallop_inf_id\'/omicsMap.csv, case clear comma
sort identifier
merge 1:1 identifier using \`scallop_inf_id\'/work/INTERVALdata_28FEB2020
format Olink_inf_QC_24m %15.0g
format Olink_inf_gwasQC_24m %15.0g
format Affymetrix_QC_bl %15.0g
format Affymetrix_gwasQC_bl %15.0g
format Affymetrix_QC_24m %15.0g
format Affymetrix_gwasQC_24m %15.0g

keep identifier Olink_inf_QC_24m Olink_inf_gwasQC_24m Affymetrix_gwasQC_bl ethnicPulse agePulse sexPulse
d,f
tab ethnicPulse
di _N
outsheet if Olink_inf_QC_24m!=. | Olink_inf_gwasQC_24m!=. using idmap.tsv,noquote replace
gen str20 ethnic=ethnicPulse
replace ethnic="EUR" if inlist(ethnicPulse,"Eng/W/Scot/NI/Brit","White Irish")==1
replace ethnic="EAS" if inlist(ethnicPulse,"Asian- Bangladeshi","Asian- Indian","Asian- Pakistani","Chinese")==1
replace ethnic="MID" if ethnicPulse=="Arab"
gen ethnic_NA=ethnic
replace ethnic_NA="NA" if inlist(ethnic,"EUR","EAS","MID")==0
gen FID=0
rename Affymetrix_gwasQC_bl IID
keep if Olink_inf_QC_24m!=. | Olink_inf_gwasQC_24m!=.
outsheet FID IID ethnic using ethnic.txt if IID!=., noquote replace
outsheet FID IID ethnic_NA using ethnic_NA.txt if IID!=., noquote replace
END

cut -f1 ${gwas_id} | grep -f - idmap.tsv | wc -l

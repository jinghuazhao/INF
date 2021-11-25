// INF1.proxy.avoutput-rsid
insheet using "INF1.proxy.avoutput-rsid", case name clear
keep SNP location FuncensGene GenerefGene
sort SNP
rename SNP Signals
save INF1.proxy, replace
// ST6
insheet using ST6., case names clear
drop v8 v9
drop if replicate=="1" | replicate =="Not found"
duplicates tag Protein Signals, gen(dup)
list Signals UniProt Protein dup
sort Protein Signals
save st6, replace
outsheet Protein Signals using ST6, noquote replace
// INF1.merge.cis.vs.trans-rsid
insheet using INF1.merge.cis.vs.trans-rsid, clear case delim(" ") names
sort prot SNP
rename prot Protein
rename SNP Signals
merge 1:1 Protein Signals using st6
keep if _merge==1
gen ccis=.
replace ccis=0 if cis=="TRUE"
replace ccis=1 if cis=="FALSE"
keep uniprot Protein Signals pprot ccis
sort Signals
merge Signals using INF1.proxy
keep if _merge==3
drop _merge pprot
gen no=_n
export excel st6.xlsx, firstrow(variables) replace
// CVD I
insheet using ST6., case names clear
keep if Source=="Folkersen, et al (2020)"
list Signals UniProt Protein Source replicate
tab Protein
// Lothian
insheet using ST6., case names clear
keep if Source=="Hillary, et al. (2019)"
list Signals UniProt Protein Source replicate
tab Protein
// AGES
insheet using ST6., case names clear
keep if Source=="Emilsson, et al. (2018)"
list Signals UniProt Protein Source replicate
tab Protein
// NSPHS
insheet using ST6., case names clear
keep if Source=="Höglund, et al (2019)"
list Signals UniProt Protein Source replicate
tab Protein


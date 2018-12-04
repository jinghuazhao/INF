# 4-12-2018 JHZ

source tryggve/analysis.ini

declare -a proteins=(ARTN CCL25 CD6 CST5 FGF.5 IFN.gamma IL.13 IL.18R1 IL.1.alpha IL.20 IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN TSLP IL.10RA IL.5 TNF)
for s in $(seq 1 ${#cohorts[@]}); do
    export protein=${cohorts[$(( $s-1 ))]}
    echo $protein
    if [ ! -d $protein ]; then mkdir $protein; fi
    R --no-save -q < tryggve/QCGWAS.R > $protein.log
done

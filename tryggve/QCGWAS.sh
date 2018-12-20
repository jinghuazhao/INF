# 20-12-2018 JHZ

source tryggve/analysis.ini
export rt=$HOME/INF

cd $rt
parallel --env rt -j1 -C' ' '
echo {}; \
export protein={}; \
cd $rt/work; \
if [ ! -d {} ]; then mkdir {}; fi; \
cd -; \
R --no-save -q < $rt/tryggve/QCGWAS.R > $rt/work/{}.log' ::: IFN.gamma IL.22.RA1 TSLP

declare -a prot=(ARTN CCL25 CD6 CST5 FGF.5 IFN.gamma IL.13 IL.18R1 IL.1.alpha IL.20 IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN IL.10RA IL.5 TNF TSLP)
function protein_array ()
{
# mkdir
  for s in $(seq 1 ${#prot[@]}); do
      export protein=${prot[$(( $s-1 ))]}
      if [ ! -d work/$protein ]; then mkdir work/$protein; fi
  done
# QCGWAS
  export protein=$1
  R --no-save -q < tryggve/QCGWAS.R > work/$protein.log
}

# 27-3-2019 JHZ

source tryggve/analysis.ini

export rt=$HOME/INF
cd $rt
parallel --env rt -j1 -C' ' '
  echo {}; \
  if [ ! -d {} ]; then mkdir $rt/work/{}; fi; \
  export protein={}; \
  R --no-save -q < $rt/tryggve/QCGWAS.R > $rt/work/{}.log
' ::: ADA CCL13 CCL25 CCL4 CCL8 CD6 CDCP1 CST5 CXCL6 FGF.5 GDNF IL.10RB IL.12B IL.15RA IL.18R1 MIP.1.alpha MMP.1 MMP.19 uPA

# 1. proteins with an alarmingly large number of signals:
# ARTN CCL25 CD6 CST5 FGF.5 IFN.gamma IL.13 IL.18R1 IL.1.alpha IL.20 IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN IL.10RA IL.5 TNF TSLP
# 2. proteins with the highest lod%:
# IFN.gamma IL.22.RA1 TSLP
# 3. proteins with a large number of signals which are likely to be genuine:
# ADA CCL13 CCL25 CCL4 CCL8 CD6 CDCP1 CST5 CXCL6 FGF.5 GDNF IL.10RB IL.12B IL.15RA IL.18R1 MIP.1.alpha MMP.1 MMP.10 uPA

function protein_array ()
# alternative implementation
{
  declare -a p=(ADA CCL13)
  for s in $(seq 1 ${#p[@]}); do
    # mkdir
      export protein=${p[$(( $s-1 ))]}
      if [ ! -d work/$p ]; then mkdir work/$p; fi
    # QCGWAS
      R --no-save -q < tryggve/QCGWAS.R > work/$p.log
  done
}

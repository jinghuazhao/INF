# 4-12-2018 JHZ

source tryggve/analysis.ini

export rt=$HOME/INF

# work on tryggve and work subdirectories

cd $rt/work
parallel --env rt -j1 -C' ' '
echo {}; \
export protein={}; \
if [ ! -d {} ]; then mkdir {}; fi; \
R --no-save -q < $rt/tryggve/QCGWAS.R > {}.log' ::: CCL25 CD6 CST5 FGF.5 IL.13 IL.18R1 IL.1.alpha IL.20 IL.20RA IL.22.RA1 IL.24 IL.2RB IL.33 LIF MCP.2 NRTN IL.10RA IL.5 TNF

## ARTN IFN.gamma TSLP

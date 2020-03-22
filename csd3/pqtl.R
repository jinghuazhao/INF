#!/rds/user/jhz22/hpc-work/bin/Rscript --vanilla

catalogue <- "pQTL"
ps <- function(rsid)
      phenoscanner::phenoscanner(snpquery=rsid, catalogue=catalogue, proxies = "EUR", pvalue = 1e-07, r2= 0.6, build=37)
allpQTL <- scan("work/INF1.merge.snp",what="")
Folkersen <- c("rs17610659", "rs5744249", "rs344562", "rs1883832", "rs574044675", "rs79287178", "rs2468187", "rs8951")
pqtl17 <- c("rs59950280", "rs5745687", "rs6475938")
suhre17 <- c("rs4815244")
sun18 <- c("rs2228145","rs704", "rs8064426", "rs113010081", "rs6827617", "rs635634", "rs579459", "rs2247769", "rs4734879", "rs9469127",
           "rs6993770", "rs705379", "rs757973", "rs13229619", "rs6921438", "rs9272226", "rs11759846", "rs10733789", "rs10822155",
           "rs7090111", "rs16840522", "rs11220462", "rs7310615", "rs3184504", "rs7137828", "rs597808", "rs653178", "rs28929474",
           "rs8178824", "rs78357146", "rs516316", "rs516246", "rs601338", "rs798893", "rs5030044", "rs3733402", "rs66530140", "rs3130510",
           "rs7763262", "rs11713634", "rs7612912", "rs73133996", "rs712048", "rs2270297", "rs16850073", "rs6073958", "rs1366949",
           "rs450373", "rs7406661", "rs34790908", "rs1133763", "rs12191307", "rs11212636", "rs12075", "rs3014874", "rs11265493",
           "rs671623", "rs3136676", "rs3859189", "rs838131", "rs7564243", "rs1260326", "rs7624160", "rs9815073", "rs1491961", "rs35728689",
           "rs4241577", "rs10076557", "rs28377109", "rs2228467", "rs62292952", "rs2032887", "rs4251805", "rs17860955", "rs12588969",
           "rs55744193", "rs471994", "rs34557412", "rs11548618", "rs2229092")
knownpQTL <- c(Folkersen,pqtl17,suhre17,sun18)
rsid <- setdiff(allpQTL,knownpQTL)
m <- length(rsid)
if (m<=100) { r <- ps(rsid) } else
{
  r1 <- ps(rsid[1:100]); lapply(r1,dim)
  r2 <- ps(rsid[101:m]); lapply(r2,dim)
  snps <- rbind(with(r1,snps),with(r2,snps))
  results <- rbind(with(r1,results),with(r2,results))
  r <- list(snps=snps,results=results)
}
lapply(r,dim)
save(r,file=paste0("work/INF1.merge.",catalogue))
options(width=500)
attach(r)
results <- within(results,{
   a1 <- ref_a1
   a2 <- ref_a2
   swap <- ref_a1 > ref_a2
   a1[swap] <- ref_a2[swap]
   a2[swap] <- ref_a1[swap]
   ref_snpid <- paste0(ref_hg19_coordinates,"_",a1,"_",a2)
})
for(d in unique(with(results,dataset)))
{
  cat(d,"\n")
  sink(paste(catalogue,d,sep="."))
  s <- subset(results,dataset==d)
  print(s[c("ref_rsid","ref_snpid","rsid","r2","p","trait")])
  sink()
}
detach(r)

## Novel pQTL 
# --- Folkerson et al. (2017) ---
# rs8951 -- MIP.1.alpha
# rs344562 -- TNFSF14
# rs1883832 -- CD40
# rs574044675 -- TRAIL
# rs79287178 -- TRANCE
# rs2468187 -- TRANCE
# --- pQTL ---
# rs59950280 -- HGF
# --- Sun et al. (2018) ---
# rs579459 -- 
# rs4734879 -- 
# rs6993770 -- CXCL5
# rs705379 -- SCF
# rs13229619 -- FGF.21
# rs9272226 -- CDCP1
# rs11759846 -- IL.1.alpha
# rs10733789 -- CCL11
# rs10822155 -- VEGF.A
# rs16840522 -- TRAIL
# rs11220462 -- uPA
# rs7310615 -- TNFB
# rs3184504 -- CXCL11, CXCL10, CD5, IL.12B, CXCL9, CD244
# rs7137828 -- MIP.1.alpha
# rs597808 -- CD6
# rs28929474 -- TRAIL
# rs8178824 -- TRAIL
# rs5030044 -- TRAIL
# rs78357146 -- IL18R1
# rs516316 -- MMP.10, CST5
# rs798893 -- SCF
# rs3733402 -- 4E.BP1
# rs66530140 -- ST1A1
# rs3130510 -- IL.12B
# rs7763262 -- CX3CL1
# rs11713634 -- TRANCE
# rs7612912 -- MCP.4
# rs73133996 -- TWEAK
# rs16850073 -- CXCL6
# rs6073958 -- SCF
# rs1366949 -- CXCL1
# rs7406661 -- uPA
# rs1133763 -- MCP.2
# rs12191307 -- CXCL9
# rs11212636 -- Flt3L
# rs12075 -- MCP.1, CCL11, MCP.2, MCP.3, CXCL6, MCP.4
# rs3014874 -- EN.RAGE
# rs11265493 -- CD244
# rs671623 -- CX3CL1
# rs3136676 -- MCP.4
# rs3859189 -- OSM
# rs838131 -- FGF.21
# rs7564243 -- uPA
# rs1260326 -- FGF.21
# rs7624160 -- Flt3L
# rs9815073 -- IL.12B
# rs1491961 -- CCL11
# rs35728689 -- MCP.1
# rs4241577 -- CXCL9
# rs10076557 -- IL.12B
# rs28377109 -- IL.10
# rs2228467 -- MCP.1, CCL11, MCP.3
# rs62292952 -- CCL19
# rs17860955 -- MMP.10
# rs12588969 -- IL.12B
# rs471994 -- MMP.1
# rs34557412 -- TNFRSF9
# rs2229092 -- TNFB

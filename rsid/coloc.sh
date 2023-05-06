#!/usr/bin/bash
cat << 'EOL' > ${INF}/work/st14.tsv
Protein	ID	Disease	P	FDR	pQTL	Cognate disease QTL	Secondary disease signal 	Disease P-val	r2 (pQTL-disease QTL)	Elimination	Reason
IL-12B	ebi-a-GCST004132	Crohn's disease	1.2E-21	3.4E-19	rs10076557	rs10046001	TRUE	6.2E-22	0.88	FALSE
CD40	ebi-a-GCST004132	Crohn's disease	2.2E-08	1.2E-06	rs1883832	rs6032664	TRUE	1.3E-07	0.99	FALSE
IL-18R1	ebi-a-GCST004132	Crohn's disease	1.9E-05	6.2E-04	rs2270297	rs11378157	FALSE	1.2E-13	0.95	FALSE
IL-18R1	ieu-a-996	Eczema	2.1E-10	1.3E-08	rs2270297	rs6419573	FALSE	2.9E-10	0.98	FALSE
IL-12B	ebi-a-GCST004131	Inflammatory bowel disease	1.5E-30	8.2E-28	rs10076557	rs10045431	TRUE	4.4E-32	0.88	FALSE
CD6	ebi-a-GCST004131	Inflammatory bowel disease	2.1E-07	1.1E-05	rs2074227	rs11230563	FALSE	2.0E-06	0.99	FALSE
CD40	ebi-a-GCST004131	Inflammatory bowel disease	1.9E-06	8.2E-05	rs1883832	rs6074022 	TRUE	1.4E-06	0.99	FALSE
CD40	ieu-b-18	Multiple sclerosis	1.2E-12	9.8E-11	rs1883832	rs6032662	FALSE	2.8E-13	0.99	FALSE
CD5	ieu-a-1112	Primary sclerosing cholangitis	8.1E-05	2.4E-03	rs674379	rs10792302	FALSE	4.9E-05	0.99	FALSE
CD40	ieu-a-833	Rheumatoid arthritis	1.4E-15	1.3E-13	rs1883832	rs4239702	FALSE	9.0E-15	0.84	FALSE
IL-12B	ebi-a-GCST004133	Ulcerative colitis	4.7E-20	6.5E-18	rs10076557	rs983825	TRUE	2.1E-20	0.99	FALSE
CXCL5	ebi-a-GCST004133	Ulcerative colitis	2.3E-06	9.0E-05	rs450373	rs425535	FALSE	2.2E-07	0.99	FALSE
EOL

Rscript -e '
  options(width=200)
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
  st14 <- read.delim(file.path(INF,"work","st14.tsv")) %>%
          left_join(pQTLdata::inf1[c("prot","target.short","gene")],by=c("Protein"="target.short"))
  res_formatted <- list()
  for(index in 1:12)
  {
    rm(d1,d2)
    Protein <- st14[index,"Protein"]
    prot <- st14[index,"prot"]
    efo <- st14[index,"ID"]
    label <- st14[index,"Disease"]
    pgwas <- read.table(file.path(INF,"mr","gsmr","prot",paste0(prot,".gz")),header=TRUE) %>%
             left_join(read.table(file.path(INF,"mr","gsmr","prot",paste0(prot,"-rsid.txt")),col.names=c("SNP","rsid")))
    d1 <- with(pgwas,list(beta=b,varbeta=se^2,N=N,MAF=if_else(freq<0.5,freq,1-freq),type="quant",snp=SNP))
    tgwas <- read.table(file.path(INF,"mr","gsmr","trait",paste0(prot,"-",efo,".gz")),header=TRUE) %>%
             mutate(freq=if_else(freq==".",0.5,as.numeric(freq))) %>%
             rename(ALT=A1,REF=A2) %>%
             left_join(select(pgwas,SNP,A1,A2)) %>%
             mutate(sw=if_else(A1==ALT,1,-1),b=sw*b) %>%
             filter(!is.na(b) & !is.na(freq))
    d2 <- with(tgwas,list(beta=b,varbeta=se^2,N=N,MAF=if_else(freq<0.5,freq,1-freq),type="quant",snp=SNP))
    coloc_res <- coloc::coloc.abf(dataset1=d1, dataset2=d2, p1=1e-4, p2=1e-4, p12=1e-5)
    res_formatted[[index]] <- data.frame(Protein,efo,label,dplyr::as_tibble(t(as.data.frame(coloc_res$summary))))
    print(res_formatted[[index]])
  }
  res <- as.data.frame(do.call(rbind,res_formatted)) %>%
         setNames(c("Protein","ID","Disease","nsnps","PP0","PP1","PP2","PP3","PP4"))
  write.table(res,file=file.path(INF,"coloc","gwas-coloc.tsv"),row.names=FALSE,quote=FALSE,sep="\t")
  library(scales)
  for(col in 5:9) res[,col] <- percent(res[,col])
'

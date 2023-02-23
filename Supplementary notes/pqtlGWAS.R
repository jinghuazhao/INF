options(width=200)

suppressMessages(library(dplyr))
suppressMessages(library(gap))
suppressMessages(library(pQTLtools))

inf1_prot <- vector()
for(i in 1:92) inf1_prot[inf1[i,"prot"]] <- mutate(inf1[i,],target.short=if_else(!is.na(alt_name),alt_name,target.short))[["target.short"]]
INF <- Sys.getenv("INF")

signs <- c(-1,0,1)
symbols <- c("-","0","+")
unicodes <- c("\u25E8","\u2605","\u21F5","\u26AB","\u25EF")
ctcols <- c("red","blue")

INF1_metal <- read.delim(file.path(INF,"work","INF1.METAL"),as.is=TRUE) %>%
              mutate(hg19_coordinates=paste0("chr",Chromosome,":",Position)) %>%
              rename(INF1_rsid=rsid, Total=N) %>%
              left_join(pQTLdata::inf1[c("prot","gene","target.short","alt_name")]) %>%
              mutate(target.short=if_else(!is.na(alt_name),alt_name,target.short),
                     pqtl_direction=sapply(Effect, function(x) {symbols[signs==sign(x)]})) %>%
              select(-alt_name)
INF1_aggr <- INF1_metal %>%
  select(Chromosome,Position,target.short,gene,hg19_coordinates,MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,log.P.,cis.trans,INF1_rsid) %>%
  group_by(Chromosome,Position,MarkerName,INF1_rsid,hg19_coordinates) %>%
  summarise(nprots=dplyr::n(),
            prots=paste(target.short,collapse=";"),
            Allele1=paste(toupper(Allele1),collapse=";"),
            Allele2=paste(toupper(Allele2),collapse=";"),
            EAF=paste(Freq1,collapse=";"),
            Effects=paste(Effect,collapse=";"),
            SEs=paste(StdErr,collapse=";"),
            directions=paste(sapply(Effect, function(x) {symbols[signs==sign(x)]}),collapse=";"),
            log10P=paste(log.P.,collapse=";"),
            cistrans=paste(cis.trans,collapse=";")) %>% data.frame()
INF1_aggr_save <- INF1_aggr
setwd(INF)
rsid <- INF1_aggr[["INF1_rsid"]]
catalogue <- "GWAS"
proxies <- "EUR"
p <- 5e-8
r2 <- 0.8
build <- 37
metal <- subset(within(INF1_metal,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
                select=-c(Chromosome,Position,Direction))
aggr <- subset(within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
               select=-c(Chromosome,Position,INF1_rsid))
imd_diseases <- read.table(file.path(INF,"ebi","efo-3.26.0","efo_0000540.csv"),col.names=c("efo","disease"),as.is=TRUE,sep=",")
efo_diseases <- read.table(file.path(INF,"ebi","efo-3.26.0","efo_diseases.csv"),col.names=c("efo","disease"),as.is=TRUE,sep=",")
iid_diseases <- read.delim(file.path(INF,"doc","immune.efos.txt"),as.is=TRUE) %>%
                mutate(id=gsub("_",":",id))

if (FALSE) {
# queries by block as in pQTLtools
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  snps_results <- with(r,right_join(snps,results))
  ps <- subset(snps_results,select=-c(hg38_coordinates,ref_hg38_coordinates,pos_hg38,ref_pos_hg38,dprime))
  sarcoidosis <- with(ps,grep("sarcoidosis",trait))
  ps[sarcoidosis,"efo"] <- "Orphanet_797"
  save(INF1_aggr,r,ps,file=file.path(INF,"ps","INF1.merge.GWAS"))
# bulk results from the Web interface
  snps <- read.delim(file.path(INF,"ps","SNP.tsv"))
  ref_rsid <- pull(snps,ref_rsid) %>% unique
  write.table(ref_rsid,file=file.path(INF,"ps","INF1_ref_rsid.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
  results <- read.delim(file.path(INF,"ps","GWAS.tsv")) # %>% mutate(chr=gsub("chr","",lapply(strsplit(results$hg19_coordinates,":"),"[",1)))
  ps <- right_join(snps,results)
# inidividual SNP queries from R
  r <- snpqueries(ref_rsid, block_size=1, waiting_time=1, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  ps <- with(r,right_join(snps,results))
# individual SNP queries from the command-line.
  snps <- read.delim(file.path(INF,"ps","SNP-single.tsv"))
  results <- read.delim(file.path(INF,"ps","GWAS-single.tsv"))
  ps <- right_join(snps,results)
  save(INF1_aggr,r,ps,file=file.path(INF,"ps","INF1.merge.GWAS"))
} else load(file.path(INF,"ps","INF1.merge.GWAS"))

ps_gcst <- mutate(ps,a1_ps=a1,a2_ps=a2)
# PMID 23128233 (heatmap has a gray cell/NA) which requires API, but gwasrapidd cames handy
# PhenoScanner list rs601338 a1/a2 (A/G)
# P=2.1e-09 is recongnised and the following gives [0.086-0.169] though treated as continuous outcome
rs601338 <- gwasrapidd::get_associations(variant_id = "rs601338")
filter(rs601338@risk_alleles,association_id=="100045013")
print(gap::ci2ms("0.086-0.169")) # b/SE=-2.12/0.17
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs601338","a1"] <- "A"
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs601338","direction"] <- "-"
# rs516246 is done similarly but there are many duplicates from GWAS Catalog so we resort to paper/ST2
# P=1e-15, CI=1.071-1.143, b/SE=0.101/0.0166
rs516246 <- gwasrapidd::get_associations(variant_id = "rs516246")
filter(rs516246@risk_alleles,association_id=="35848328")
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs516246","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs516246","direction"] <- "+"
# PhenoScanner has multiple entries but OpenGWAS is ready for check REF/ALT (C/G) ES:SE:LP=0.00143911:0.000261552:7.42556
# https://gwas.mrcieu.ac.uk/files/ukb-a-100/ukb-a-100.vcf.gz
ps_gcst[ps_gcst$pmid=="UKBB" & ps_gcst$rsid=="rs7310615" & ps_gcst$p=="3.754e-08","a1"] <- "G"
ps_gcst[ps_gcst$pmid=="UKBB" & ps_gcst$rsid=="rs7310615" & ps_gcst$p=="3.754e-08","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="26151821" & ps_gcst$rsid=="rs3184504","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="26151821" & ps_gcst$rsid=="rs3184504","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="26502338" & ps_gcst$rsid=="rs597808","a1"] <- "A"
ps_gcst[ps_gcst$pmid=="26502338" & ps_gcst$rsid=="rs597808","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="27992413" & ps_gcst$rsid=="rs3184504","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="27992413" & ps_gcst$rsid=="rs3184504","direction"] <- "+"
# PhenoScanner and GCST agree
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs11230563","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="23128233" & ps_gcst$rsid=="rs11230563","direction"] <- "+"
# PhenoScanner and GCST disagree on risk allele -- the latter has "A" (in A/G) in line with the paper (on - strand)
ps_gcst[ps_gcst$pmid=="22961000" & ps_gcst$rsid=="rs3184504" & ps_gcst$efo=="EFO_1001486","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="22961000" & ps_gcst$rsid=="rs3184504" & ps_gcst$efo=="EFO_1001486","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="23603763" & ps_gcst$rsid=="rs3184504" & ps_gcst$efo=="EFO_0004268","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="23603763" & ps_gcst$rsid=="rs3184504" & ps_gcst$efo=="EFO_0004268","direction"] <- "+"
# PhenoScanner provides a guess than GCST
ps_gcst[ps_gcst$pmid=="27182965" & ps_gcst$rsid=="rs635634","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="27182965" & ps_gcst$rsid=="rs635634","direction"] <- "+"
# PhenoScanner/GCST definitions (G/T) of risk allele are different
# but the paper indicates T as minor allele and using PLINK therefore GCST is correct
ps_gcst[ps_gcst$pmid=="22672568" & ps_gcst$rsid=="rs495828","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="22672568" & ps_gcst$rsid=="rs495828","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="24262325" & ps_gcst$rsid=="rs579459","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="24262325" & ps_gcst$rsid=="rs579459","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="25939597" & ps_gcst$rsid=="rs7725218","a1"] <- "G"
ps_gcst[ps_gcst$pmid=="25939597" & ps_gcst$rsid=="rs7725218","direction"] <- "+"
# According to rs3184504 risk/other (T/C),
# PhenoScanner is correct to indicate rs653178 risk allele to be "C" but rs3184504 should have risk allele "T"
# since LDhap (https://ldlink.nci.nih.gov/?tab=ldhap) indicates rs3184504/rs653178 TC=0.5278, r2=0.9449
# However, GWAS catalogue appears to be wrong about beta/se since it treats Tonsillectomy OR as log(OR)
#'         where 0.035-0.067 has beta/se=0.0507/0.0081 rather than -3.03/0.166
# The GWAS Catalog is correctly in line with Nat Comm paper.
# PhenoScanner has C allele to be +
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs635634","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs635634","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs3184504","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs3184504","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs653178","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs653178","direction"] <- "-"
# PhenoScanner has multiple conflicting entries
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs681343","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs681343","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs516316","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="28928442" & ps_gcst$rsid=="rs516316","direction"] <- "+"
# T1D
# Read from PLoS Genet paper according to P value, fixing a1 by GWAS catalog and adding direction
# GWAS Catalog indicates OR=1.3(4) but misses beta/se
ps_gcst[ps_gcst$pmid=="21829393" & ps_gcst$rsid=="rs3184504","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="21829393" & ps_gcst$rsid=="rs3184504","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="22493691" & ps_gcst$rsid=="rs3184504","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="22493691" & ps_gcst$rsid=="rs3184504","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="23603761" & ps_gcst$rsid=="rs7137828","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="23603761" & ps_gcst$rsid=="rs7137828","direction"] <- "+"
# See correction below
ps_gcst[ps_gcst$pmid=="24390342" & ps_gcst$rsid=="rs1950897","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="24390342" & ps_gcst$rsid=="rs1950897","direction"] <- "+"
# correction needed on PhenoScanner, https://www.nature.com/articles/ng.789
#       snp      rsid a1   beta direction                     trait     pmid
# rs1950897  rs911263  C 0.2546         + Primary biliary cirrhosis 21399635
ps_gcst[ps_gcst$pmid=="21399635" & ps_gcst$rsid=="rs911263","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="21399635" & ps_gcst$rsid=="rs911263","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="28029757" & ps_gcst$rsid=="rs28929474","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="28029757" & ps_gcst$rsid=="rs28929474","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="25802187" & ps_gcst$rsid=="rs1883832","a1"] <- "T"
ps_gcst[ps_gcst$pmid=="25802187" & ps_gcst$rsid=="rs1883832","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="22232737" & ps_gcst$rsid=="rs3784099","a1"] <- "A"
ps_gcst[ps_gcst$pmid=="22232737" & ps_gcst$rsid=="rs3784099","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="26192919" & ps_gcst$rsid=="rs516246" & ps_gcst$efo=="EFO_0000384","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="26192919" & ps_gcst$rsid=="rs516246" & ps_gcst$efo=="EFO_0000384","direction"] <- "+"
ps_gcst[ps_gcst$pmid=="26192919" & ps_gcst$rsid=="rs516246" & ps_gcst$efo=="EFO_0003767","a1"] <- "C"
ps_gcst[ps_gcst$pmid=="26192919" & ps_gcst$rsid=="rs516246" & ps_gcst$efo=="EFO_0003767","direction"] <- "+"

ps_gsub <- ps_gcst %>%
           mutate(trait=gsub("Coronary artery disease age 50|Coronary artery disease males","Cardiovascular diseases",trait)) %>%
           mutate(trait=if_else(unit=="-"&grepl("Arthritis rheumatoid|Rheumatoid arthritis",trait),"rheumatoid arthritis",trait)) %>%
           mutate(trait=if_else(unit=="-"&grepl("Diabetes mellitus type 1|Type 1 diabetes",trait),"Type I diabetes",trait)) %>%
           mutate(trait=if_else(unit=="-"&grepl("Inflammatory bowel disease",trait),"inflammatory bowel disease",trait)) %>%
           mutate(trait=gsub("including oligoarticular and rheumatoid factor negative polyarticular JIA","",trait)) %>%
           mutate(trait=gsub("Blood clot in the leg|Self-reported deep venous thrombosis","deep vein thrombosis",trait)) %>%
           mutate(trait=gsub("Blood clot in the lung","pulmonary embolism",trait)) %>%
           mutate(trait=gsub("Celiac disease","Coeliac disease",trait)) %>%
           mutate(trait=gsub("Coronary heart disease","Coronary artery disease",trait)) %>%
           mutate(trait=gsub("Crohns","Crohn's",trait)) %>%
           mutate(trait=gsub("Advanced age related macular degeneration","Age-related macular degeneration",trait)) %>%
           mutate(trait=gsub("Renal overload gout|Renal underexcretion gout","Gout",trait)) %>%
           mutate(trait=gsub("Hayfever, allergic rhinitis or eczema|Allergic disease|Allergic disease asthma hay fever or eczema","allergy",trait)) %>%
           mutate(trait=gsub("High grade serous ovarian cancer|Invasive ovarian cancer","ovarian cancer",trait)) %>%
           mutate(trait=gsub("Serous invasive ovarian cancer|Serous boarderline ovarian cancer","ovarian cancer",trait)) %>%
           mutate(trait=gsub("Other rheumatoid arthritis","rheumatoid arthritis",trait)) %>%
           mutate(trait=gsub("Self-reported cholelithiasis or gall stones","cholelithiasis",trait)) %>%
           mutate(trait=gsub("erythematosis","erythematosus",trait)) %>%
           mutate(trait=gsub(" or myxoedema","",trait)) %>%
           mutate(trait=gsub(" or thyrotoxicosis","",trait)) %>%
           mutate(trait=gsub("Type 2 diabetes","Type II diabetes",trait)) %>%
           mutate(trait=gsub(" or endometrial","",trait)) %>%
           mutate(trait=gsub(" [+] or [-] dvt","",trait)) %>%
           mutate(trait=gsub("Illnesses of siblings: ","",trait)) %>%
           mutate(trait=gsub("Selective IgA deficiency","Selective IgA deficiency disease",trait)) %>%
           mutate(trait=gsub("Doctor diagnosed |heart attack or |Self-reported |Illnesses of father: |Illnesses of mother: ","",trait)) %>%
           mutate(trait=gsub("high blood pressure","hypertension",trait)) %>%
           mutate(trait=gsub("malabsorption or |Low grade and borderline serous |Mouth or teeth dental problems: ","",trait)) %>%
           mutate(trait=gsub("Unspecified |Vascular or heart problems diagnosed by doctor: ","",trait)) %>%
           mutate(trait=gsub("Acute myocardial infarction|Myocardial infarction|myocardial infarction","cardiovascular diseases",trait)) %>%
           mutate(trait=gsub("Coronary artery disease|heart attack|heart disease","cardiovascular diseases",trait)) %>%
           mutate(trait=gsub(" or sle| or large artery stroke| or ischemic stroke","", trait)) %>%
           mutate(trait=gsub("\\b(^[a-z])","\\U\\1",trait,perl=TRUE))
ps_filter <- ps_gsub %>%
             filter(!grepl("None of the above",trait)) %>%
             filter(!(rsid=="rs635634" & snp=="rs579459" & grepl("Cardiovascular diseases|High cholesterol|Hypertension|Type II diabetes",trait))) %>%
             filter(!(snp=="rs635634" & grepl("Deep vein thrombosis|Pulmonary embolism",trait) & proxy=="1")) %>%
             filter(!(snp=="rs579459" & grepl("Pulmonary embolism|Phlebitis and thrombophlebitis|Deep vein thrombosis|Ovarian cancer|Phlebitis and thrombophlebitis",trait) & proxy=="1")) %>%
             filter(!(snp=="rs579459" & (grepl("Cardiovascular diseases|Hypertension",trait) & proxy=="0" | direction=="NA"))) %>%
             filter(!(rsid=="rs579459" & snp=="rs635634" & grepl("High cholesterol|Ovarian cancer|Phlebitis and thrombophlebitis",trait))) %>%
             filter(!(rsid=="rs597808" & efo=="EFO_0004705;EFO_1001055" & proxy=="1")) %>%
             filter(!(rsid=="rs597808" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003777",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs516246" & pmid!="26192919")) %>%
             filter(!(rsid=="rs66530140" & grepl("Deep vein thrombosis|Pulmonary embolism",trait) & proxy=="1")) %>%
             filter(!(rsid=="rs653178" & grepl("Hypertension",trait) & proxy=="1")) %>%
             filter(!(rsid=="rs653178" & grepl("EFO_0000378|EFO_0000612|EFO_0001060",efo) & direction=="NA")) %>%
             filter(!(rsid=="rs653178" & grepl("EFO_0001060|EFO_0003956;EFO_0005854;EFO_0000274|EFO_0004705;EFO_1001055",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs653178" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003767|EFO_0003777",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs7137828" & grepl("EFO_0003785|EFO_0004705",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs7137828" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003777",efo) & proxy=="1")) %>%
             filter(!(rsid %in% c("rs516246","rs516316") & grepl("ring up phlegm",trait))) %>%
             filter(!(rsid %in% c("rs516246","rs516316") & grepl("EFO_0000384|EFO_0004799|EFO_0008111",efo) & proxy=="1")) %>%
             filter(!(efo%in%c("EFO_0004211","EFO_0004530") & direction=="NA")) %>%
             filter(!(rsid=="rs7310615" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003777",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs3184504" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003777",efo) & proxy=="0" & direction=="NA")) %>%
             filter(!(rsid=="rs3184504" & grepl("EFO_0000612|EFO_0000378|EFO_1000883|EFO_0001645|EFO_0003777",efo) & proxy=="1")) %>%
             filter(!(rsid=="rs3184504" & efo=="EFO_0004268" & proxy=="1")) %>%
             filter(!(rsid=="rs3184504" & grepl("EFO_0004705",efo) & proxy=="1")) %>%
             filter(!(rsid %in% c("rs3184504","rs597808","rs7137828","rs7310615") & grepl("EFO_0003956|EFO_0005854|EFO_0000274",efo) & proxy=="1")) %>%
             filter(!(rsid %in% c("rs3184504","rs653178") & efo=="EFO_0000537" & direction=="NA")) %>%
             filter(!(rsid %in% c("rs516246","rs516316","rs597808","rs3184504","rs7137828","rs7310615") & efo=="EFO_0000537" & proxy=="1")) %>%
             filter(!(rsid %in% c("rs597808","rs3184504","rs7137828","rs7310615") & efo=="EFO_0004325" & proxy=="1")) %>%
             filter(!(pmid=="22961000" & efo=="EFO_1001486" & proxy=="1")) %>%
             filter(!(pmid=="UKBB" & rsid=="rs7310615" & efo=="EFO_0004705;EFO_1001055" & proxy=="1")) %>%
             filter(!(pmid=="21829393" & proxy==1 & trait=="Type 1 diabetes")) %>%
             filter(!(pmid=="26151821" & proxy==1 & trait=="Colorectal cancer")) %>%
             filter(!(pmid=="27992413" & proxy==1 & trait=="Primary sclerosing cholangitis")) %>%
             filter(!(pmid=="26192919" & rsid=="rs516246" & proxy==1 & trait=="Inflammatory bowel disease")) %>%
             filter(!(pmid=="26192919" & rsid=="rs601338" & trait=="Inflammatory bowel disease")) %>%
             filter(!(pmid=="UKBB" & efo=="EFO_0000676" & proxy==1)) %>%
             filter(!(rsid=="rs601338" & grepl("Crohn's disease|Cholelithiasis|Bring up phlegm|High cholesterol|Hypertension",trait) & proxy=="1")) %>%
             filter(!(rsid=="rs601338" & efo=="EFO_0008111" & proxy=="1")) %>%
             filter(!(rsid=="rs601338" & pmid=="23128233")) %>%
             filter(!(rsid=="rs601338" & pmid=="26192919" & proxy=="1")) %>%
             filter(!(rsid=="rs7137828" & (direction=="NA" | proxy=="1" & grepl("Juvenile",trait)))) %>%
             filter(!(snp=="rs601338" & pmid=="22482804")) %>%
             filter(!(pmid=="23143596")) %>%
             filter(!(pmid=="24390342" & efo=="EFO_0000685" & is.na(beta))) %>%
             filter(!(pmid=="21907864" & is.na(beta))) %>%
             filter(!(pmid=="25646370")) %>%
             filter(!(pmid=="20383146" & direction=="NA")) %>%
             filter(!(pmid=="20081858")) %>%
             filter(!(pmid=="27790247")) %>%
             filter(!(pmid=="22399527" & direction=="NA")) %>%
             filter(!(pmid=="21386085")) %>%
             filter(!(pmid=="20453842" & direction=="NA")) %>%
             filter(!(pmid=="19430483"|pmid=="27117709"|pmid=="27197191")) %>%
             filter(!(pmid=="25305756"|pmid=="20167578"|pmid=="27997041")) %>%
             filter(!(pmid=="27182965"|pmid=="27618447"|pmid=="21383967"|pmid=="22057235")) %>%
             filter(!(pmid=="22561518"|pmid=="21383967")) %>%
             filter(!(pmid=="18794853"|pmid=="26752265"|pmid=="19430480"|pmid=="28067908")) %>%
             filter(!(pmid=="21383967"|pmid=="26621817"|pmid==""|pmid==""|pmid=="")) %>%
             filter(!(pmid=="20453842" & direction=="NA")) %>%
             filter(!(pmid=="26691988" & direction=="NA")) %>%
             filter(!(pmid=="22399527" & direction=="NA")) %>%
             filter(!(pmid=="21378990" & direction=="NA")) %>%
             filter(!(pmid=="22672568" & direction=="NA")) %>%
             filter(!(pmid=="22434691" & direction=="NA")) %>%
             filter(!(pmid=="21399635" & direction=="NA")) %>%
             filter(!(pmid=="24390342" & direction=="NA")) %>%
             filter(!(pmid=="22493691" & (direction=="NA"|proxy=="1"))) %>%
             filter(!(pmid=="26192919" & direction=="NA")) %>%
             filter(!(pmid=="20190752" & direction=="NA")) %>%
             filter(!(pmid=="27723758" & direction=="NA")) %>%
             filter(!grepl("INVT|IVNT|SDS|Z-score|bpm|crease|g/l|kg|lu|mg|ml|mmHg|mol|years|ug|unit|%",unit)) %>%
             filter(!(unit=="-"&(pmid=="UKBB"|grepl("Cholesterol ldl|Intercellular adhesion molecule 1",trait)))) %>%
             filter(!(unit=="-"&grepl("Protein quantitative trait loci|Receptors interleukin 6|Monocyte chemoattractant protein 1",trait))) %>%
             filter(!(unit=="-"&grepl("Blood proteins|Chemokine ccl2|Monocyte chemoattractant protein 1|Uric acid",trait))) %>%
             filter(!(unit=="-"&grepl("Alkaline phosphatase|selectin|Lipid metabolism|Vascular endothelial growth factor a",trait))) %>%
             filter(!(unit=="-"&grepl("Cholesterol|Cholesterol hdl|Vitamin B12|Cystatin c|Interleukin 18|Phospholipids|Metabolism",trait))) %>%
             filter(!(unit=="-"&grepl("Inflammatory biomarkers|Paraoxonase activity|Waist circumference and related phenotypes",trait))) %>%
             filter(!(unit=="-"&grepl("Age related disease endophenotypes|Age related diseases mortality and associated endophenotypes",trait))) %>%
             filter(!(unit=="-"&grepl("Glucose transporter type 2|Glucose tolerance test|Type 1 diabetes autoantibodies",trait))) %>%
             filter(!(unit=="-"&grepl("Angiotensin converting enzyme inhibitors|Celiac disease or Rheumatoid arthritis",trait))) %>%
             filter(!(unit=="-"&grepl("Coronary artery disease or ischemic stroke|Coronary artery disease or large artery stroke",trait))) %>%
             filter(!(unit=="-"&grepl("kidney stone or ureter stone/bladder stone",trait))) %>%
             filter(!trait=="Self-reported diabetes") %>%
             filter(!trait=="Self-reported eczema or dermatitis") %>%
             filter(!trait=="Illnesses of father: chronic bronchitis or emphysema") %>%
             filter(!grepl("Thyroid peroxidase antibody positivity|Smoking status: previous",trait)) %>%
             filter(!grepl("Started insulin within one year diagnosis of diabetes|No blood clot",trait)) %>%
             filter(!grepl("Medication for pain relief|Pain type experienced in last month",trait)) %>%
             filter(!grepl("Qualifications: college or university degree",trait)) %>%
             filter(!grepl("Medication for cholesterol, blood pressure or diabetes:",dataset)) %>%
             filter(!grepl("Astle|GIANT|GLGC|MAGIC",dataset)) %>%
             filter(!grepl("Illnesses of siblings: none|Treatment|Taking other prescription medications",trait)) %>%
             filter(!grepl("Body Mass index|Blood pressure|Fibrinogen|Fasting|Gamma glutamyltransferase",trait)) %>%
             filter(!grepl("Eosinophils|Gamma glutamyltransferase|mass index",trait)) %>%
             filter(!grepl("Lipid metabolism phenotypes|Triglycerides|C reactive protein",trait)) %>%
             filter(!grepl("Diabetes diagnosed by doctor|Types of physical activity|Cause of death: alcoholic hepatic failure",trait)) %>%
             filter(!grepl("Medication for cholesterol, blood pressure or diabetes|No treatment with medication for cholesterol",trait)) %>%
             filter(!grepl("Neovascularization No treatment with medication|Autism spectrum disorder or schizophrenia",trait)) %>%
             filter(!grepl("Qualifications: none",trait)) %>%
             filter(!grepl("Vascular or heart problems diagnosed by doctor: none of the above",trait)) %>%
             filter(!grepl("count|density|education|intake|levels|weight",trait)) %>%
             filter(!(dataset=="GRASP" & !grepl("Celiac disease|Coronary artery disease|Hypertension|Hypothyroidism|JIA|Myocardial infarction|Primary biliary cirrhosis|Primary sclerosing cholangitis|Type 1 diabetes|Generalized vitiligo|Rheumatoid arthritis and celiac disease",trait))) %>%
             filter(!(pmid=="21980299" & efo=="EFO_0001359")) %>%
             rename(disease=trait)
#          mutate(trait=gsub("\\b(^[A-Z])","\\L\\1",trait,perl=TRUE))
ps_mutate <- ps_filter

overlap <- function(dat,f1,f2)
# Now on individual proteins
{
  mat <- select(dat,prot,target.short,gene,hgnc,snp,MarkerName,INF1_rsid,cis.trans,Effect,StdErr,pqtl_direction,Allele1,Allele2,
                    chr,rsid,a1,a2,a1_ps,a2_ps,efo,ref_rsid,ref_a1,ref_a2,proxy,r2,
                    HLA,beta,se,p,direction,disease,n_cases,n_controls,unit,ancestry,pmid,study) %>%
         mutate(a2=if_else(a1==a1_ps,a2_ps,a1_ps)) %>%
         mutate(prefix=if_else(HLA==1,paste0(gene,"-",INF1_rsid,"-",cis.trans,"*"),paste0(gene,"-",INF1_rsid,"-",cis.trans)),
                rsidProt=paste0(prefix," (",hgnc,")"), Trait=gsub("\\b(^[a-z])","\\U\\1",disease,perl=TRUE),
         #      Effect=round(Effect,3), StdErr=round(StdErr,3),
                r2=round(as.numeric(r2),3),
         #      beta=round(as.numeric(beta),3), se=round(as.numeric(se),3), p=format(as.numeric(p),digits=3,scientific=TRUE),
                Allele1=toupper(Allele1), Allele2=toupper(Allele2), a1=toupper(a1), a2=toupper(a2),
                ref_rsid_a1_a2=paste0(ref_rsid,":",ref_a1,":",ref_a2),
                snp_rsid_chr=paste(snp,rsid,chr,sep="_")) %>%
         #      left_join(haps) %>%
         mutate(hap=paste0(Allele1,a1),pah=paste0(Allele2,a2),
                h11=paste0(Allele1,a1),h12=paste0(Allele1,a2),h21=paste0(Allele2,a1),h22=paste0(Allele2,a2),
                switch=case_when(proxy==0 & Allele1==a1 ~ "0", proxy==0 & Allele1!=a1 ~ "1",
                                 proxy==1 & hap==h11 & pah %in% c(h12,h21,h22) ~ "0", proxy==1 & hap==h12 & pah %in% c(h11,h21,h22) ~ "1",
                                 proxy==1 & hap==h21 & pah %in% c(h11,h12,h22) ~ "1", proxy==1 & hap==h22 & pah %in% c(h11,h12,h21) ~ "0",
                                 TRUE ~ "0"),
                direction=case_when(switch=="1" & direction=="-" ~ "+", switch=="1" & direction=="+" ~ "-", TRUE ~ direction),
                pqtl_trait_direction=paste0(pqtl_direction,direction),
                trait_direction=case_when(pqtl_trait_direction=="++" ~ "1",  pqtl_trait_direction=="+-" ~ "-1",
                                          pqtl_trait_direction=="-+" ~ "-1", pqtl_trait_direction=="--" ~ "1",
                                          pqtl_trait_direction=="-NA" ~ "NA",
                                          TRUE ~ as.character(direction))) %>%
       # filter(direction%in%c("-","+")) %>%
         select(-c(snp_rsid_chr,hap,pah,h11,h12,h21,h22)) %>%
         filter(!grepl("^NA",rsidProt))
  combined <- group_by(mat,hgnc,rsidProt,Trait,desc(n_cases)) %>%
              summarize(ref_rsids=paste(ref_rsid_a1_a2,collapse=";"),
                        proxies=paste(proxy,collapse=";"),
                        directions=paste(trait_direction,collapse=";"),
                        units=paste(unit,collapse=";"),
                        studies=paste(study,collapse=";"),
                        PMIDs=paste(pmid,collapse=";"),
                        diseases=paste(disease,collapse=";"),
                        efos=paste(efo,collapse="+")
                       ) %>% data.frame()
  efo_Traits <- with(combined,unique(Trait))
  rsid_Prots <- with(combined,unique(rsidProt))
  rxc <- matrix(0, length(efo_Traits), length(rsid_Prots), dimnames=list(efo_Traits,rsid_Prots))
  dn <- matrix("", length(efo_Traits), length(rsid_Prots), dimnames=list(efo_Traits,rsid_Prots))
  for(rn in rownames(rxc)) for(cn in colnames(rxc)) {
     cnrn <- subset(combined,Trait==rn & rsidProt==cn)
     if(nrow(cnrn)==0) next
     val <- unlist(strsplit(cnrn[["directions"]],";"))
     tab <- table(val)
     sym <- sapply(unlist(strsplit(cnrn[["directions"]],";"))[1],function(x){signs[symbols==x]})
     if (length(tab)==1) rxc[rn,cn] <- as.numeric(val[1]) else dn[rn,cn] <- unicodes[1]
  }
  # all beta's are NAs when unit=="-"
  subset(mat[c("study","pmid","unit","beta","pqtl_direction","direction")],unit=="-")
  # all studies with risk difference were UKBB
  subset(mat[c("study","pmid","unit","beta","n_cases","n_controls","pqtl_direction","direction")],unit=="risk diff")
  write.table(select(mat,-prot,-MarkerName,-prefix,-rsidProt,-INF1_rsid,-a1_ps,-a2_ps,-Effect,-StdErr,rsid,-beta,-se,-p,
                         -n_cases,-n_controls,-pqtl_trait_direction,-trait_direction,-Trait) %>%
              rename(Protein=target.short,Target_gene=gene,Nearest_gene=hgnc,Proxy=proxy,EFO=efo,Disease=disease,PMID=pmid,Study=study),
              file=file.path(INF,"ps",f1),row.names=FALSE,quote=FALSE,sep=",")
  write.table(select(combined,-desc.n_cases.),file=file.path(INF,"ps",f2),row.names=FALSE,quote=FALSE,sep=",")
  list(rxc=rxc,dn=dn)
}

SF <- function(rxc, dn, f="SF-pQTL-IMD-GWAS.png", ch=21, cw=21, h=16, w=17, ylab="Immune-mediated outcomes")
{
  print(rownames(rxc))
  print(dim(rxc))
  library(grid)
  library(pheatmap)
  col <- colorRampPalette(c("#4287f5","#ffffff","#e32222"))(3)
  png(file.path(INF,"ps",f),res=300,width=w,height=h,units="in")
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  p <- pheatmap(rxc, legend=FALSE, angle_col="315", border_color="black", color=col, cellheight=ch, cellwidth=cw,
                 display_numbers=dn, number_color = "brown",
                 cluster_rows=TRUE, cluster_cols=TRUE, fontsize=16)
  ccols <- if_else(grepl("-cis",p$gtable$grobs[[4]]$label),1,2)
  print(cbind(ccols,ctcols[ccols],p$tree_col$labels,p$gtable$grobs[[4]]$label))
  p$tree_col$labels <- gsub("^[0-9]*-|-cis|-trans","",p$tree_col$labels)
  p$gtable$grobs[[4]]$label <- gsub("^[0-9]*-|-cis|-trans","",p$gtable$grobs[[4]]$label)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  colnames(rxc) <- gsub("^[0-9]*-|-cis|-trans","",colnames(rxc))
  p <- pheatmap(rxc, legend=FALSE, angle_col="315", border_color="black", color=col, cellheight=ch, cellwidth=cw,
                 display_numbers=dn, number_color = "brown",
                 cluster_rows=TRUE, cluster_cols=TRUE, fontsize=24)
  p$gtable$grobs[[4]]$gp=gpar(col=ctcols[ccols])
  setHook("grid.newpage", NULL, "replace")
  grid.draw(p)
  grid.text("Protein-pQTL (Nearest gene)", y=0.01, gp=gpar(fontsize=28))
  grid.text(ylab, x=-0.01, rot=90, gp=gpar(fontsize=28))
  dev.off()
}

# GWAS diseases
long <- merge(metal,ps_mutate,by="hg19_coordinates",all.y=TRUE)
dat <- long
f1 <- "ST-pQTL-disease-overlap.csv"
f2 <- "ST-pQTL-disease-overlap-combined.csv"
rxc_gwas <- overlap(dat,f1,f2)
with(rxc_gwas,SF(rxc,dn,f="SF-pQTL-disease-overlap.png",ch=35,cw=35,h=42,w=52,ylab="GWAS diseases"))

# All EFOs for IMD but somehow smaller number of rows
imd_list <- imd_diseases[["efo"]]
imd_list <- iid_diseases[["id"]]
sel <- sapply(gsub("_",":",long[["efo"]]),function(x) {
             long_set <- unlist(strsplit(x,";"))
             set_int <- intersect(long_set,imd_list)
             paste0(imd_list[imd_list%in%set_int],collapse="")!=""
             })
dat <- filter(long,sel)
f1 <- "ST-pQTL-IMD-overlap.csv"
f2 <- "ST-pQTL-IMD-overlap-combined.csv"
rxc_imd2 <- overlap(dat,f1,f2)
with(rxc_imd2,SF(rxc,dn,f="SF-pQTL-IMD-overlap.png",ch=35,cw=35,h=20,w=36))

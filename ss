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
imd_diseases <- read.table(file.path(INF,"ebi","efo-3.26.0","efo_0000540.csv"),col.names=c("efo","disease"),as.is=TRUE,sep=",") %>%
                mutate(efo=gsub(":", "_", efo))
efo_diseases <- read.table(file.path(INF,"ebi","efo-3.26.0","efo_diseases.csv"),col.names=c("efo","disease"),as.is=TRUE,sep=",") %>%
                mutate(efo=gsub(":", "_", efo))
iid_diseases <- read.delim(file.path(INF,"doc","immune.efos.txt"),as.is=TRUE)

if (FALSE) {
  r <- snpqueries(rsid, catalogue=catalogue, proxies=proxies, p=p, r2=r2, build=build)
  lapply(r,dim)
  snps_results <- with(r,right_join(snps,results))
  ps <- subset(snps_results,select=-c(hg38_coordinates,ref_hg38_coordinates,pos_hg38,ref_pos_hg38,dprime))
  sarcoidosis <- with(ps,grep("sarcoidosis",trait))
  ps[sarcoidosis,"efo"] <- "Orphanet_797"
  INF1_aggr <- INF1_aggr_save
  save(INF1_aggr,r,ps,file=file.path(INF,"work","INF1.merge.GWAS"))
} else {
  load(file.path(INF,"work","INF1.merge.GWAS"))
  INF1_aggr <- INF1_aggr_save
}
ps <- ps %>%
      mutate(a1_ps=a1,a2_ps=a2)
metal <- subset(within(INF1_metal,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
                select=-c(Chromosome,Position,Direction))
aggr <- subset(within(INF1_aggr,{HLA <- as.numeric(Chromosome==6 & Position >= 25392021 & Position <= 33392022)}),
               select=-c(Chromosome,Position,INF1_rsid))
immune_infection_efo <- with(iid_diseases,gsub(":","_",id))
# crude processing
ps_na_direction <- subset(ps,direction=="NA") %>%
                   filter(!grepl("FVII|HDL|IFT172|IgG|LDL|Albumin|Fasting|NFATotal|Plasma|Serum|Triglycerides|Vitamin",trait)) %>%
                   filter(!grepl("levels|plasma|reactive",trait))
write.table(unique(sort(ps_na_direction[["trait"]])),file=file.path(INF,"work","ps_na_direction.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)
# manually entered and can be modified lexicographically
na_selected <- c(
"Advanced age related macular degeneration",
"Age related disease",
"Allergy",
"Alopecia areata",
"Antineutrophil cytoplasmic antibody associated vasculitis",
"Arthritis rheumatoid",
"Asthma",
"Autism spectrum disorder or schizophrenia",
"Breast cancer",
"Breast cancer",
"Cancer pleiotropy",
"Celiac disease",
"Childhood ear infection",
"Chronic hepatitis B infection",
"Coronary artery disease",
"Depression",
"Diabetes mellitus type 1",
"Extreme obesity with early age of onset",
"Generalized vitiligo",
"Glaucoma primary open angle",
"Gout",
"Hypertension",
"Hypothyroidism",
"Idiopathic membranous nephropathy",
"IgA nephropathy",
"Inflammatory bowel disease",
"Insulin resistance",
"Ischemic stroke",
"Juvenile idiopathic arthritis including oligoarticular and rheumatoid factor negative polyarticular JIA",
"Kidney diseases",
"Liver cirrhosis biliary",
"Metabolic syndrome x",
"Myocardial infarction",
"Primary biliary cirrhosis",
"Primary sclerosing cholangitis",
"Prostate cancer",
"Rheumatoid arthritis",
"Rheumatoid arthritis cyclic citrullinated peptide CCP positive",
"Schizophrenia",
"Selective IgA deficiency",
"Systemic lupus erythematosus SLE",
"Tonsillectomy",
"Type 1 diabetes",
"Type 2 diabetes",
"Venous thrombo",
"Vitiligo")

ps_na_grep <- paste(na_selected,collapse="|")
ps_na_disease <- filter(ps_na_direction,grepl(ps_na_grep,trait))
ps_na_disease[c("snp","rsid","proxy","r2","trait","pmid","study")]
ps_na_pmid <- unique(ps_na_disease$pmid) %>%
              paste(collapse=" ")
write.table(ps_na_disease,file=file.path(INF,"work","ps_na_disease.csv"),quote=FALSE,row.names=FALSE,sep=",")

load(file.path(INF,"work","GCST.rda"))
ps_na_gcst <- function(rsid,PMID)
{
  publications <- select(GCST_studies@publications, study_id, pubmed_id) %>%
                  filter(pubmed_id %in% PMID)
  studies <- select(GCST_studies@studies, study_id, reported_trait) %>%
             filter(study_id %in% publications$study_id)
  sources <- left_join(publications,studies)
  risk_alleles <- select(GCST@risk_alleles, association_id, variant_id, risk_allele, risk_frequency) %>%
                  filter(variant_id==rsid)
  associations <- select(GCST@associations, association_id, range, beta_unit, beta_direction, beta_number, standard_error, pvalue)
  id <- filter(assoc_study,association_id %in% associations$association_id) %>% distinct()
  r <- left_join(id,publications) %>% filter(!is.na(pubmed_id)) %>% distinct()
  risk_allels <- risk_alleles %>% filter(association_id %in% r$association_id) %>% distinct()
  associations <- associations %>% filter(association_id %in% r$association_id) %>% distinct()
# get_variants(variant_id=rsid,pubmed_id=PMID)
  left_join(r,associations) %>% left_join(risk_alleles) %>% data.frame() %>% select(-association_id,-study_id)
}

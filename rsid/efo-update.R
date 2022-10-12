suppressMessages(library(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(gwasrapidd))

# Input
options(width=200)
HPC_WORK <- Sys.getenv("HPC_WORK")
INF <- Sys.getenv("INF")
rt <- file.path(INF,"mr","gsmr")

xlsx <- "https://jhz22.user.srcf.net/INF/latest/efo-update.xlsx"
efo_update <- read.xlsx(xlsx,sheet=1) %>%
              filter(!is.na(Disease)) %>%
              mutate(opengwasid=unlist(lapply(strsplit(Source,"/"),"[",5)))
efo_info <- ieugwasr::gwasinfo(pull(efo_update,opengwasid)) %>%
            mutate(trait=gsub("\\b(^[a-z])","\\U\\1",substr(trait,1,30),perl=TRUE)) %>%
            select(id,trait,ncase,ncontrol,nsnp,population,pmid,author,year)
knitr::kable(efo_info)
efo_add <- read.xlsx(xlsx,sheet=2)
knitr::kable(efo_add)

csv <- read.delim(file.path(INF,"OpenGWAS","finngen_endpoints.tsv"))
sel.trait <- c("D3_SARCOIDOSIS","M13_ANKYLOSPON","M13_SJOGREN","AB1_HIV","AB1_VIRAL_MENINGITIS")
sel.var <- c("phenotype","phenocode","number.of.cases","number.of.controls")
finngen <- subset(csv[sel.var],phenocode %in%sel.trait) %>%
           rename(id=phenocode,trait=phenotype,ncase=number.of.cases,ncontrol=number.of.controls) %>%
           mutate(id=paste0("finngen_R7_",id))
knitr::kable(subset(csv,phenocode %in%sel.trait)[sel.var] %>% arrange(phenotype))

efo_update <- bind_rows(efo_info,efo_add,finngen)
write.table(efo_update,file=file.path(INF,"OpenGWAS","efo-update.txt"),quote=FALSE,row.names=FALSE,sep="\t")

check_efo <- function(efo_id,out,trait.name=FALSE)
{
  if (trait.name) id <- get_studies(efo_trait=efo_id)
  else id <- get_studies(trait_to_study(efo_id) %>% pull(study_id))
  write.table(slot(id,"publications") %>%
              data.frame() %>%
              filter(!grepl("UK Biobank|Tobacco",title) & publication_date > "2017-01-09") %>%
              select(-author_orcid),file=out,sep="\t",row.names=FALSE,quote=FALSE)
}

gwas_catalog_check <- function()
{
  xlsx <- "https://jhz22.user.srcf.net/INF/latest/efo.xlsx"
  efo <- read.xlsx(xlsx,sheet=1,startRow=2) %>%
         filter(grepl("ukb|bbj|finn",MRBASEID) & is.na(Replacement)) %>%
         select(-c(uri,Zhengetal,Replacement, Cases, Controls)) %>%
         arrange(MRBASEID)

  gs <- get_studies(trait_to_study(pull(efo,id)) %>% pull(study_id))
  gs@ancestries %>% data.frame()

  check_efo("EFO_0003767","ibd.tsv") # no change
  check_efo("Orphanet_797","sarcoidosis.tsv") # no change
  check_efo("EFO_0004237","graves.tsv") # only BBJ
  check_efo("EFO_0000764","hiv.tsv") # GCST90096801/2 only with hits
  check_efo("EFO_0000699","sjorgren.tsv") # GCST004062 only with partial information
  check_efo("EFO_0000584","meningitis.tsv") # unavailable
  check_efo("EFO_0004208","vitiligo.tsv") # GCST010676 only with hits
  check_efo("EFO_0003106","pneumonia.tsv") # GCST90134363 10k hits, https://research.23andme.com/dataset-access/
  check_efo("Orphanet_3389","tuberculosis.tsv") # GCST004922/3, only with hits
  check_efo("EFO_0005854","allergic-rhinitis.tsv") # only UKB
  check_efo("EFO_0003785","allergy.tsv") # unavailable
  check_efo("EFO_0004705","hypothyroidism.tsv") # GCST006898 only with hits or UKB
  check_efo("EFO_0003103","urinary-tract-infection.tsv") # only UKB
# ---
  check_efo("EFO_0000270","asthma.tsv") # unavailable since replaced with http://purl.obolibrary.org/obo/MONDO_0004979
  check_efo("MONDO_0004979","asthma.tsv")
  check_efo("EFO_0004591","asthma.tsv") # childhood asthma replaced with http://purl.obolibrary.org/obo/MONDO_0005405
  check_efo("MONDO_0005405","asthma.tsv")
  check_efo("asthma","asthma.tsv",trait.name=TRUE)
  # https://ftp.ebi.ac.uk/pub/databases/gwassummary_statistics/GCST90131001-GCST90132000/ (AA)
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005212/ (TAGC)
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010043/ (UKB+TAGC)
  # based on OpenGWAS for GCST90014325
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014325/
  check_efo("EFO_0002609","jia.tsv") # BBJ paper, GCST90010715
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90010001-GCST90011000/GCST90010715/
  check_efo("EFO_1001486","pbc.tsv") # EUR/AS considerably larger, GCST90061440[1|2] but only one entry availabhle
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90061001-GCST90062000/GCST90061440/
  check_efo("EFO_0001060","celiac.tsv") # GCST005319 [21|22], GCST009874, only with hits
  check_efo("EFO_0001359","t1d.tsv") # unavailable
  # http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014023/
  check_efo("EFO_0000270","asthma.tsv") # unavailable
  check_efo("EFO_0004194","iga.tsv") # GCST011296, summary data unavailable; # GCST011781 only with hits
  check_efo("EFO_0001068","malaria.tsv") # GCST009514[5] summary statistics unavailable
  check_efo("EFO_0000764","hiv.tsv") # GCST90096801[2] summary statistics unavailable
  check_efo("EFO_1001231","uveitis.tsv") # GCST006126[7] only with hits, BBJ
  check_efo("MONDO_0100096","covid.tsv") # GCST90134599,GCST90134598,GCST90134602,GCST90134601,GCST90134596,GCST90134597,GCST90134600
  # COVID-19 hospitalisation
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90134001-GCST90135000/GCST90134602/
  check_efo("EFO_0000676","psoriasis.tsv") # GCST90019016
  # https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90019001-GCST90020000/GCST90019016/
  check_efo("EFO_0003898","ankyspondy.tsv") # no useful information
  check_efo("EFO_0004268","scleroschol.tsv") # no useful information
  check_efo("EFO_0003779","thyroid.tsv") # BBJ
  check_efo("EFO_0000584","meningitis.tsv") # no data

  immune_infection <- read.delim(file.path(INF,"doc","immune.efos.txt"),as.is=TRUE) %>%
                      mutate(id=gsub(":","_",id)) %>%
                      filter(! id %in% c("EFO_0003840","EFO_0003885"))
}

# https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/

efo_2022_10_05 <- function()
{
  xlsx <- "https://jhz22.user.srcf.net/INF/latest/GSMR_datasets.xlsx"
  efo <- read.xlsx(xlsx,sheet=1,startRow=2,colNames=TRUE,skipEmptyRows=TRUE)
  efo_left <- efo[1:5] %>%
              rename(Source1=Source) %>%
              mutate(N.cases=as.numeric(N.cases),N.controls=as.numeric(N.controls),Total.N=as.numeric(Total.N))
  efo_right <- efo[6:10] %>% rename(Source2=Source) %>% mutate(N.cases=as.numeric(N.cases),N.controls=as.numeric(N.controls),Total.N=as.numeric(Total.N))
  efo_old <- filter(cbind(efo_left,efo_right[5]),!grepl("http",Source2)) %>%
             select(-Source2) %>% rename(Source=Source1)
  efo_new <- filter(efo_right,grepl("http",Source2)) %>%
             rename(Source=Source2)
  efo_update <- bind_rows(efo_old,efo_new) %>%
                filter(!is.na(Disease))
  missed <- with(efo_update, grepl("Crohn",Disease))
  efo_update[missed,c("N.cases","N.controls","Source")] <- c(12194,280722,"https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST004132/")
  ra <- with(efo_update, grepl("Rheumatoid",Disease))
  efo_update[ra,c("N.cases","N.controls","Source")] <- c(19234,61565,"https://gwas.mrcieu.ac.uk/datasets/ieu-a-833")
  efo_update <- mutate(efo_update,opengwasid=unlist(lapply(strsplit(Source,"/"),"[",5)))
  efo_update <- arrange(efo_update,opengwasid)
  sel.trait <- c("JUVEN_ARTHR","D3_SARCOIDOSIS","M13_SLE")
  finngen <- with(efo_update,grepl("finn-b",opengwasid))
  efo_update[finngen,c("N.cases","N.controls","Total.N")]
  finngen_r7_N <- subset(csv,phenocode %in%sel.trait)[sel.var] %>%
                  arrange(phenocode) %>%
                  select(number.of.cases, number.of.controls) %>%
                  mutate(all=number.of.cases+number.of.controls)
  efo_update[finngen,c("N.cases","N.controls","Total.N")] <- finngen_r7_N
  write.table(efo_update,file="efo-update.txt",quote=FALSE,row.names=FALSE,sep="\t")
  knitr::kable(efo_update)
  # Output
  xlsx <- file.path("efo-update.xlsx")
  wb <- createWorkbook(xlsx)
  hs <- createStyle(textDecoration="BOLD", fontColour="#FFFFFF", fontSize=12, fontName="Arial Narrow", fgFill="#4F80BD")
  for (sheet in c("efo_update"))
  {
     addWorksheet(wb,sheet,zoom=150)
     writeData(wb,sheet,sheet,xy=c(1,1),headerStyle=createStyle(textDecoration="BOLD",
               fontColour="#FFFFFF", fontSize=14, fontName="Arial Narrow", fgFill="#4F80BD"))
     body <- get(sheet)
     writeDataTable(wb, sheet, body, xy=c(1,2), headerStyle=hs, firstColumn=TRUE, tableStyle="TableStyleMedium2")
     freezePane(wb, sheet, firstCol=TRUE, firstActiveRow=3)
     width_vec <- apply(body, 2, function(x) max(nchar(as.character(x))+2, na.rm=TRUE))
   # width_vec_header <- nchar(colnames(body))+2
     setColWidths(wb, sheet, cols = 1:ncol(body), widths = width_vec)
     writeData(wb, sheet, tail(body,1), xy=c(1, nrow(body)+2), colNames=FALSE, borders="rows", borderStyle="thick")
  }
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}

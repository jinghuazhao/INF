#!/usr/bin/bash

function coloc()
{
cat << 'EOL' > ${INF}/work/st14.tsv
Protein	ID	Disease	P	FDR	pQTL	Cognate disease QTL	Secondary disease signal 	Disease P-val	r2 (pQTL-disease QTL)	Elimination	Reason
IL-12B	ebi-a-GCST004132	Crohns disease	1.2E-21	3.4E-19	rs10076557	rs10046001	TRUE	6.2E-22	0.88	FALSE
CD40	ebi-a-GCST004132	Crohns disease	2.2E-08	1.2E-06	rs1883832	rs6032664	TRUE	1.3E-07	0.99	FALSE
IL-18R1	ebi-a-GCST004132	Crohns disease	1.9E-05	6.2E-04	rs2270297	rs11378157	FALSE	1.2E-13	0.95	FALSE
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
          left_join(pQTLdata::inf1[c("prot","target.short","gene")],
                    by=c("Protein"="target.short"))
  res_formatted <- list()
  for(index in 1:12)
  {
    rm(d1,d2)
    Protein <- st14[index,"Protein"]
    prot <- st14[index,"prot"]
    efo <- st14[index,"ID"]
    label <- st14[index,"Disease"]
    pgwas <- read.table(file.path(INF,"mr","gsmr","prot",paste0(prot,".gz")),header=TRUE) %>%
             left_join(read.table(file.path(INF,"mr","gsmr","prot",paste0(prot,"-rsid.txt")),
                                  col.names=c("SNP","rsid")))
    d1 <- with(pgwas,list(beta=b,varbeta=se^2,N=N,MAF=if_else(freq<0.5,freq,1-freq),
               type="quant",snp=SNP))
    tgwas <- read.table(file.path(INF,"mr","gsmr","trait",paste0(prot,"-",efo,".gz")),
                        header=TRUE) %>%
             mutate(freq=if_else(freq==".",0.5,as.numeric(freq))) %>%
             rename(ALT=A1,REF=A2) %>%
             left_join(select(pgwas,SNP,A1,A2)) %>%
             mutate(sw=if_else(A1==ALT,1,-1),b=sw*b) %>%
             filter(!is.na(b) & !is.na(freq))
    tgwas <- read.table(file.path(INF,"mr","gsmr","trait",paste0(prot,"-",efo,".gz")),
                        header=TRUE) %>%
             mutate(freq=if_else(freq==".",0.5,as.numeric(freq))) %>%
             rename(ALT=A1,REF=A2) %>%
             left_join(select(pgwas,SNP,A1,A2)) %>%
             mutate(sw=if_else(A1==ALT,1,-1),b=sw*b) %>%
             filter(!is.na(b) & !is.na(freq))
    d2 <- with(tgwas,list(beta=b,varbeta=se^2,N=N,MAF=if_else(freq<0.5,freq,1-freq),
               type="quant",snp=SNP))
    coloc_res <- coloc::coloc.abf(dataset1=d1, dataset2=d2, p1=1e-4, p2=1e-4, p12=1e-5)
    res_formatted[[index]] <- data.frame(Protein,efo,label,
                                         dplyr::as_tibble(t(as.data.frame(coloc_res$summary))))
    print(res_formatted[[index]])
  }
  res <- as.data.frame(do.call(rbind,res_formatted)) %>%
         setNames(c("Protein","ID","Disease","nsnps","PP0","PP1","PP2","PP3","PP4"))
  write.table(res,file=file.path(INF,"coloc","gwas-coloc.tsv"),
              row.names=FALSE,quote=FALSE,sep="\t")
  library(scales)
  for(col in 5:9) res[,col] <- percent(res[,col])
'
}

function cis_lst()
{
  export M=250000
  Rscript -e '
    suppressMessages(library(dplyr))
    INF <- Sys.getenv("INF")
    M <- Sys.getenv("M") %>% as.integer
    st14 <- read.delim(file.path(INF,"work","st14.tsv")) %>%
            select(Protein,ID,Disease,pQTL)
    INF1_METAL <- read.delim(file.path(INF,"work","INF1.METAL")) %>%
                  left_join(select(pQTLdata::inf1,prot,target.short,gene,chr,start,end)) %>%
                  left_join(st14,by=c("target.short"="Protein","rsid"="pQTL")) %>%
                  filter(cis.trans=="cis",!is.na(Disease)) %>%
                  mutate(start=if_else(start-M<0,0,start-M),end=end+M,
                         region=paste0(chr,":",start,"-",end)) %>%
                  select(uniprot,prot,gene,region,chr,rsid,target.short,ID,Disease)
    write.table(INF1_METAL,file=file.path(INF,"coloc","cis.lst"),
                row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  '
}

function lz()
{
  module load python/2.7
  export dir=${INF}/coloc
  cut -f5 --complement ${dir}/cis.lst | \
  awk -vFS="\t" -vOFS="\t" '{split($4,a,":|-");print $1,$2,$3,a[1],a[2],a[3],$5,$6,$7,$8,$9}' | \
  parallel -j12 --env dir -C '\t' '
  export prot_efo={1}-{2}-{3}-{9}
  export gene_efo_pqtl={3}-{9}_{7}
# GWAS
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${INF}/mr/gsmr/trait/{2}-{9}-rsid.txt \
            --delim space title="{9}-{10}-{7}" \
            --markercol SNP --pvalcol p --chr {4} --start {5} --end {6} --cache None \
            --no-date --plotonly --prefix=GWAS-{3}-{9} --rundir ${dir} --refsnp {7}
# SCALLOP/INF
  cat <(echo -e "snpid rsid chr pos a1 a2 mlog10p") \
      <(tabix ${INF}/METAL/{2}-1.tbl.gz {4}:{5}-{6} | \
        awk "
        {
          split(\$3,a,\"_\")
          print a[1],\$1,\$2,-\$12,\$3,toupper(\$4),toupper(\$5)
        }" | \
        sort -k1,1 | \
        join -12 -21 <(grep chr{4} ${INF}/work/snp_pos) - | \
        awk -vOFS="\t" "{print \$6, \$2, \$3, \$4, \$7, \$8, \$5}") | \
  gzip -f > ${dir}/INF-${prot_efo}.tsv.gz
  (
    echo -e "chr\tpos\trsid\tmlog10P"
    gunzip -c ${dir}/INF-${prot_efo}.tsv.gz | \
    awk -v OFS="\t" "NR>1 {print \$3,\$4,\$2,\$7}" | \
    sort -k1,1n -k2,2n
  ) > ${dir}/INF-${prot_efo}.lz
  locuszoom --source 1000G_Nov2014 --build hg19 --pop EUR --metal ${dir}/INF-${prot_efo}.lz \
            --delim tab title="SCALLOP: {3}-{7}" \
            --markercol rsid --pvalcol mlog10P --no-transform --chr {4} --start {5} --end {6} --cache None \
            --no-date --plotonly --prefix=INF-{3}-{9} --rundir ${dir} --refsnp {7}
  rm ${dir}/INF-${prot_efo}.lz
  for src in GWAS INF
  do
     pdftopng -f 1 -l 1 -r 300 ${dir}/${src}-${gene_efo_pqtl}.pdf ${dir}/${src}-${gene_efo_pqtl}
     mv ${dir}/${src}-${gene_efo_pqtl}-000001.png ${dir}/${src}-${gene_efo_pqtl}.png
  done
    convert -append ${dir}/GWAS-${gene_efo_pqtl}.png ${dir}/INF-${gene_efo_pqtl}.png \
            -resize x500 -density 300 ${dir}/combine-${gene_efo_pqtl}.png
    convert ${dir}/combine-${gene_efo_pqtl}.png -quality 0 ${dir}/combine-${gene_efo_pqtl}.jp2
    convert ${dir}/combine-${gene_efo_pqtl}.jp2 ${dir}/combine-${gene_efo_pqtl}.pdf
    rm -f ${dir}/combine-${gene_efo_pqtl}.png
    rm -f ${dir}/combine-${gene_efo_pqtl}.jp2
    rm ${dir}/GWAS-${gene_efo_pqtl}.png
    rm ${dir}/INF-${gene_efo_pqtl}.png
  '
  qpdf --empty --pages $(ls ${dir}/combine*.pdf) -- ${dir}/protein-disease-lz.pdf
}

function run_PWCoCo()
{
export dir=${INF}/coloc
cat <<'EOL'> ${dir}/pwcoco.sb
#!/usr/bin/bash

#SBATCH --job-name=_pwcoco
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio

#SBATCH --export ALL
#SBATCH --output=DIR/_pwcoco.o
#SBATCH --error=DIR/_pwcoco.e

module load ceuadmin/PWCoCo/1.0
export dir=DIR

cut -f5 --complement ${dir}/cis.lst | \
awk -vFS="\t" -vOFS="\t" '{split($4,a,":|-");print $1,$2,$3,a[1],a[2],a[3],$5,$6,$7,$8,$9}' | \
parallel -j1 --env dir -C '\t' '
  export gene_efo_pqtl={3}-{9}_{7}
  export ref=~/rds/public_databases/1000G/ALL.chr{4}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  if [ ! -f ${dir}/${gene_efo_pqtl}.freq ]; then
     bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/EUR_AF\n" -r {4}:{5}-{6} ${ref} | \
     awk "{if(\$4<\$5) snpid=\"chr\"\$1\":\"\$2\"_\"\$4\"_\"\$5;
                  else snpid=\"chr\"\$1\":\"\$2\"_\"\$5\"_\"\$4;
           print snpid,\$3,\$4,\$5,\$6}" | \
     sort -k1,1 > ${dir}/${gene_efo_pqtl}.freq
  fi
  cat <(gunzip -c ${INF}/mr/gsmr/trait/{2}-{9}.gz | head -1) \
      <(gunzip -c ${INF}/mr/gsmr/trait/{2}-{9}.gz | sed "1d" | \
        join - ${dir}/${gene_efo_pqtl}.freq | \
        awk "NF>8" | \
        awk "{if(\$4==\".\") {if(\$2==\$10) \$4=1-\$12;else \$4=\$12}; print}") | \
        cut -d" " -f1-8 \
      > ${dir}/${gene_efo_pqtl}-pwcoco.sst1
  gunzip -c ${INF}/mr/gsmr/prot/{2}.gz > ${dir}/${gene_efo_pqtl}-pwcoco.sst2
  pwcoco --bfile ${INF}/INTERVAL/cardio/INTERVAL \
         --chr {4} \
         --sum_stats1 ${dir}/${gene_efo_pqtl}-pwcoco.sst1 \
         --sum_stats2 ${dir}/${gene_efo_pqtl}-pwcoco.sst2 \
         --p_cutoff1 1e-6 --p_cutoff2 5e-8 \
         --log ${dir}/${gene_efo_pqtl}-pwcoco \
         --out ${dir}/${gene_efo_pqtl}-pwcoco --out-cond
'
EOL

sed -i "s|DIR|${dir}|" ${dir}/pwcoco.sb
sbatch --wait ${dir}/pwcoco.sb

Rscript -e '
  suppressMessages(library(dplyr))
  st14 <- read.delim("~/INF/work/st14.tsv") %>%
          select(Protein,ID,Disease,pQTL)
  coloc <- data.frame()
  for (f in dir("~/INF/coloc/",pattern="-pwcoco.coloc"))
  {
     ids=unlist(strsplit(f,"-|_"))
     d <- read.delim(file.path("~/INF/coloc/",f)) %>%
          mutate(Dataset1=paste(ids[2],ids[3],ids[4],sep="-"),Dataset2=ids[1],pQTL=ids[5]) %>%
          rename(ID=Dataset1,Gene=Dataset2) %>%
          select(ID,Gene,pQTL,SNP1,SNP2,H0,H1,H2,H3,H4,log_abf_all)
     coloc <- bind_rows(coloc,d)
  }
  coloc <- left_join(st14,coloc)
  coloc.names <- names(coloc)
  H.names <- coloc.names[grepl("^H",coloc.names)]
  coloc[c(H.names,"log_abf_all")] <- round(coloc[c(H.names,"log_abf_all")],digits=2)
  write.table(coloc,file="~/INF/coloc/coloc-all.txt",row.names=FALSE,quote=FALSE,sep="\t")
  write.table(subset(coloc,H4>=0.8),file="~/INF/coloc/coloc.txt",row.names=FALSE,quote=FALSE,sep="\t")
'
}

run_PWCoCo

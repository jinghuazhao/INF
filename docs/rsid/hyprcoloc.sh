#!/usr/bin/bash

function rs12075()
{
export chr=1
export start=159175353
export end=159525679
export M=1e6
export r=$(awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} 'BEGIN{print chr":"start-M"-"end-M}')

if [ ! -d ${INF}/hyprcoloc ]; then mkdir ${INF}/hyprcoloc; fi
# Monoocyte count
# VARIANT ID_dbSNP49      CHR     BP      REF     ALT     EFFECT_INT      SE_INT  MLOG10P_INT     ALT_FREQ_INT    INFO_INT
export trait=mono
gunzip -c ~/rds/results/public/gwas/blood_cell_traits/astle_2016/raw_results/blood_cell_traits/gzipped_interval/${trait}.tsv.gz | \
awk -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} -vOFS="\t" '
{
  if(NR==1) print "chr","pos","snpid","a1","a2","af","b","se","p","n";
  else if($3==chr) print $3,$4,$2,$6,$5,$10,$7,$8,$9,173480
}' | \
bgzip -f > ${trait}.tsv.gz
module load jdk-8u141-b15-gcc-5.4.0-p4aaopt
module load gatk
module load python/3.7
cd ${HPC_WORK}/gwas2vcf
source env/bin/activate
if [ ! -f ${INF}/${trait}.vcf.gz ]; then
   python main.py --out ${INF}/${trait}.vcf.gz --data ${INF}/${trait}.tsv.gz \
                  --id ${trait} --ref human_g1k_v37.fasta --json ${INF}/rsid/gwas2vcf.json
   tabix -f ${INF}/hyprcoloc/${trait}.vcf.gz
fi
deactivate
cd -

cd ${INF}/hyprcoloc
bcftools merge -r ${r} -O z -o MCP-rs12075.vcf.gz \
         METAL/gwas2vcf/MCP.1.vcf.gz METAL/gwas2vcf/MCP.2.vcf.gz METAL/gwas2vcf/MCP.3.vcf.gz METAL/gwas2vcf/MCP.4.vcf.gz ${trait}.vcf.gz
bcftools query -f "%CHROM:%POS\_%REF\_%ALT[\t%ES][\t%SE]\n" hyprcoloc/MCP-rs12075.vcf.gz > MCP-rs12075.txt

R --no-save -q <<END
  hyprcoloc_test <- function()
  {
  # Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)
    betas <- hyprcoloc::test.betas
    head(betas)
    ses <- hyprcoloc::test.ses
    head(ses)
  # Trait names and SNP IDs
    traits <- paste0("T", 1:10)
    rsid <- rownames(betas)
  # Colocalisation analyses
    hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
  }
  require(hyprcoloc)
  hyprcoloc_test()
  n.traits <- 5
  rs12075 <- within(read.table("hyprcoloc/MCP-rs12075.txt",col.names=c("rsid",paste0("v",1:(2*n.traits)))),
             {
                 v1 <- as.numeric(v1); v2 <- as.numeric(v2)
                 v3 <- as.numeric(v3); v4 <- as.numeric(v4)
                 v5 <- as.numeric(v5); v6 <- as.numeric(v6)
                 v7 <- as.numeric(v7); v8 <- as.numeric(v8)
                 v9 <- as.numeric(v9); v10 <- as.numeric(v10)
             })
  na <- apply(is.na(rs12075),1,any)
  rs12075 <- subset(rs12075,!na)
  rsid <- with(rs12075,rsid)
  betas <- as.matrix(rs12075[paste0("v",1:n.traits)])
  ses <- as.matrix(rs12075[paste0("v",(n.traits+1):(2*n.traits))])
  colnames(betas) <- colnames(ses) <- traits <- paste0("T",1:n.traits)
  rownames(betas) <- rownames(ses) <- rsid
  hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
  d <- read.table("rs12075-beta.gassoc",col.names=c("snpid","marker","chr","pos","MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count"))
  markers <- d[c("marker","chr","pos")]
  betas <- as.matrix(d[c("MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count")])
  rownames(betas) <- with(d,marker)
  d <- read.table("rs12075-se.gassoc",col.names=c("snpid","marker","chr","pos","MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count"))
  ses <- as.matrix(d[c("MCP.1","MCP.2","MCP.3","MCP.4","Monocyte_count")])
  rownames(ses) <- with(d,marker)
  hyprcoloc(betas, ses, trait.names=colnames(ses), snp.id=with(markers,marker))
END
cd -

#Results:
#  iteration                                     traits posterior_prob
#1         1 MCP.1, MCP.2, MCP.3, MCP.4, Monocyte_count          0.978
#  regional_prob candidate_snp posterior_explained_by_snp dropped_trait
#1        0.9964       rs12075                          1            NA

zcat Whole_Blood.allpairs.txt.gz | cut -f1 | awk 'NR>1{print substr($1,1,15)}' | sort | uniq > work/wb
R --no-save -q <<END
   wb <- scan("work/wb", what="")
   g <- within(subset(grex::grex(wb),!is.na(uniprot_id)),{uniprot=trimws(uniprot_id)})
   m <- merge(gap::inf1,g,by="uniprot",all.x=TRUE)
END
}

function LTBR()
{
export chr=12
export pos=6514963
export gene=LTBR
export rsid=rs2364485
export flank_kb=1000
export b1=6300000
export b2=6700000
export bracket=${b1}-${b2}

}

#!/usr/bin/bash

function gwasglue()
{
cat << 'EOL' > ${INF}/coloc/gwasglue.sb
#!/usr/bin/bash

#SBATCH --job-name=_pwcoco
#SBATCH --mem=28800
#SBATCH --time=12:00:00

#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio

#SBATCH --export ALL
#SBATCH --output=DIR/_gwasglue.o
#SBATCH --error=DIR/_gwasglue.e

Rscript -e '

run_pwcoco <- function(id1,id2,region,
                       path_to_pwcoco="/usr/local/Cluster-Apps/ceuadmin/PWCoCo/1.0/build/pwcoco",
                       type1="quant",type2="cc",
                       workdir=Sys.getenv("HPC_WORK"))

{
  cat("id1=",id1,"id2=",id2,"\n")
  INF <- Sys.getenv("INF")
  attach(region)
  bfile <- file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr))
# region
  chrompos <- paste0(chr,":",start-M,"-",end+M)
  detach(region)
  suppressMessages(library(gwasglue))
  res <- pwcoco(id1,id2,bfile,chrompos,path_to_pwcoco,type1=type1,type2=type2,workdir=workdir)
  print(res)
}

{
  M <- 100000

# CD5 - Primary sclerosing cholangitis
  region <- list(chr=11,start=60869867,end=60895324,M=M)
  run_pwcoco("ebi-a-GCST90000452","ieu-a-1112",region)
# CD40 - Multiple sclerosis
  region <- list(chr=20,start=44746911,end=44758502,M=M)
  run_pwcoco("ebi-a-GCST90011997","ieu-b-18",region)
}

'
EOL

sed -i "s|DIR|${INF}/coloc|" ${INF}/coloc/gwasglue.sb
sbatch --wait ${INF}/coloc/gwasglue.sb
}

gwasglue

# --bfile - specifies the location of the reference dataset, normally from Plink, in the bed/bim/fam formats. Each of the bed/bim/fam files should have the same name and in the same directory.
# --sum_stats1 - first file or folder containing summary statistics. Please see below "Input" section.
# --sum_stats2 - second file or folder containing summary statistics. Please see below "Input" section.
# --log - specifies log name, default is "pwcoco_log.txt" and will save in the same folder from where the program is run.
# --out - prefix for the result files, default is "pwcoco_out".
# --p_cutoff - P value cutoff for SNPs to be selected by the stepwise selection process, default is 5e-8. Alternatively, the flags --p_cutoff1 and --p_cutoff2 may be used to specify dataset-specific P value cutoffs, relative to the order of the data given in the --sum_stats flags.
# --chr - when reading the reference files, the program will limit the analysis to those SNPs on this chromosome.
# --top_snp - maximum number of SNPs that may be selected by the stepwise selection process, default is 1e10, i.e. a lot.
# --ld_window - distance (in kb) that, when exceeded, is assumed for SNPs to be in total LE, default is 1e7.
# --collinear - threshold that, when exceeded, determines if SNPs are collinear, default is 0.9.
# --maf - filters SNPs from the reference dataset according to this threshold, default is 0.1.
# --freq_threshold - SNPs in the phenotype datasets which differ by more than this amount in the reference dataset will be excluded, default is 0.2.
# --init_h4 - PWCoCo will run an initial colocalisation on the unconditioned dataset. If the H4 for this analysis reaches this threshold, the program will terminate early. Default is 80 (i.e. 80%). Set to 0 if you would like the program to always continue regardless of the initial colocalisation result.
# --out_cond - would you like for the conditioned data to be saved as text files as well? Just including this flag will work (no extra argument following this flag is necessary).
# --coloc_pp - specify the three prior probability Ps: the next **three** arguments must be the P values, default is 1e-4, 1e-4 and 1e-5.
# --n1 - also --n2, specify the sample size (see also next flag) for the corresponding summary statistics.
# --n1_case - also --n2_case, specify the number of cases for the corresponding summary statistics.
# --threads - sets number of threads available for OpenMP multi-threaded functions, default is 8.
# --verbose - if this flag is given, PWCoCo will output files which can be used for debugging purposes. These files include SNPs which did not match the allele frequency given in the reference data and included SNPs within the analysis. Also sets `--out_cond` flag. (No extra argument following this flag is necessary).


function init()
{
  export chr=12
  export pos=6514963
  export ichip=~/rds/results/public/gwas/multiple_sclerosis/ImmunoChip_Results/Immunochip_FinalResults_LimitedDiscovery.txt
# CHR BPHG18 BPHG19 ImmunochipID Risk_Allele Ref_Allele Risk_Allele_Freq N P OR SE Q I Region
# 1 1108138 1118275 vh_1_1108138 G A 0.9586 11 7.62E-01 1.012 0.0404 0.085 39.51 none
# 1 1110294 1120431 vh_1_1110294 A G 0.0438 11 5.09E-01 1.026 0.0385 0.178 28.07 none
  awk -vchr=${chr} -vpos=${pos} -vM=100000 'NR==1||($1==chr && $3>=pos-M && $3<pos+M){
      if(NR==1) print "SNP","A1","A2","b","se","p","N";
      if($5<$6) snpid="chr"$1":"$3"_"$5"_"$6;
      else snpid="chr"$1":"$3"_"$6"_"$5
      if(NR>1) print snpid,$5,$6,$7,log($10),$11,$9,10000
  }' ${ichip} > ${INF}/MS/ichip.ma

  export v3=~/rds/results/public/gwas/multiple_sclerosis/discovery_metav3.0.meta.gz
# CHR BP SNP A1 A2 N P OR
# 1 11154 chr1:11154 C A 4 0.7911 0.9818
# 1 11565 chr1:11565 G T 5 0.8735 0.9924

  gunzip -c ${v3} | \
  awk -vchr=${chr} -vpos=${pos} -vM=100000 'NR==1||($1==chr && $2>=pos-M && $2<pos+M){
      if(NR==1) print "SNP","A1","A2","P","OR";
      if($4<$5) snpid="chr"$1":"$2"_"$4"_"$5;
      else snpid="chr"$1":"$2"_"$5"_"$4
      if(NR>1) print snpid,$4,$5,$7,$8
  }' | \
  Rscript -e '
     v3 <- within(read.table("stdin",header=TRUE),{b=log(OR);se=b/qnorm(1-P/2);EAF=0.1;N=10000});
     INF <- Sys.getenv("INF")
     write.table(v3[c("SNP","A1","A2","EAF","b","se","P","N")],file=file.path(INF,"MS","v3.ma"),quote=FALSE,row.names=FALSE)
  '

  gunzip -c ${INF}/METAL/TNFB-1.tbl.gz | \
  awk -vchr=${chr} -vpos=${pos} -vM=100000 'NR==1||($1==chr && $2>=pos-M && $2<pos+M)' | \
  cut -f3-6,10-12,18 | \
  tr '\t' ' ' | \
  awk '{if(NR>1) {$2=toupper($2);$3=toupper($3);$7=10^$7}};1' > ${INF}/MS/LTA.ma

  module load ceuadmin/PWCoCo/1.0
  pwcoco --bfile ${INF}/INTERVAL/cardio/INTERVAL \
         --chr 12 \
         --sum_stats1 ${INF}/MS/ichip.ma \
         --sum_stats2 ${INF}/MS/LTA.ma \
         --p_cutoff1 1e-6 --p_cutoff2 5e-8 \
         --log ${INF}/MS/ichip-LTA-pwcoco \
         --out ${INF}/MS/ichip-LTA-pwcoco --out-cond

  pwcoco --bfile ${INF}/INTERVAL/cardio/INTERVAL \
         --chr 12 \
         --sum_stats1 ${INF}/MS/v3.ma \
         --sum_stats2 ${INF}/MS/LTA.ma \
         --p_cutoff1 1e-6 --p_cutoff2 5e-8 \
         --log ${INF}/MS/v3-LTA-pwcoco \
         --out ${INF}/MS/v3-LTA-pwcoco
}

function deprecated()
{
awk '{print $7,$8,$9,$10,$5,$6,$11,$12}' ${INF}/MS/TNFB-pQTL-rs2364485.dat > ${INF}/MS/LTA.ma

pwcoco --bfile ${INF}/INTERVAL/cardio/INTERVAL \
       --chr 12 \
       --sum_stats1 ${INF}/MS/EUR-ieu-a-1025.ma \
       --sum_stats2 ${INF}/MS/LTA.ma \
       --p_cutoff1 1e-6 --p_cutoff2 5e-8 \
       --log ieu-a-1025-LTA-pwcoco \
       --out ieu-a-1025-LTA-pwcoco
}

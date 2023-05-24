#!/usr/bin/bash

#SBATCH --job-name=_eQTL
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --array=4,6,20,22,34,50,52
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_coloc_%A_%a.err
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_coloc_%A_%a.out
#SBATCH --export ALL

read prot gene chr start37 end37 snpid ensGene group start end < <(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]' ${INF}/coloc/cis.dat)
export prot=${prot}
export gene=${gene}
export chr=${chr}
export start37=${start37}
export end37=${end37}
export ensGene=${ensGene}
export region37=${chr}:${start37}-${end37}
export region=${chr}:${start}-${end}

function cis()
{
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(rtracklayer))
    options(width=200)
    M <- 1e6
    cvt <- "~/INF/work/INF1.merge.cis.vs.trans"
    csv <- read.csv(cvt) %>%
           filter(cis.trans=="cis") %>%
           select(prot,p.gene,p.chr,p.start,p.end,SNP) %>%
           mutate(p.start=if_else(p.start-M<0,0,p.start-M),p.end=p.end+M) %>%
           left_join(select(pQTLdata::inf1,prot,ensembl_gene_id))
    INF <- Sys.getenv("INF")
    gr <- with(csv,GenomicRanges::GRanges(seqnames=p.chr,IRanges::IRanges(p.start,p.end,names=SNP))) %>% unique()
    path <- system.file(package="pQTLtools", "eQTL-Catalogue", "hg19ToHg38.over.chain")
    ch <- import.chain(path)
    seqlevelsStyle(gr) <- "UCSC"
    gr38 <- liftOver(gr,ch) %>%
            as_tibble() %>%
            group_by(group,group_name) %>%
            summarize(start=min(start),end=max(end))
    cis <- left_join(csv,gr38,by=c("SNP"="group_name")) %>%
           data.frame
    write.table(cis,file="~/INF/coloc/cis.dat",row.names=FALSE,col.names=FALSE,quote=FALSE)
  '
}

function sumstats38()
{
  if [ ! -d ${INF}/coloc/sumstats ]; then mkdir -p ${INF}/coloc/sumstats; fi
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(rtracklayer))
    options(width=200)
    prot <- Sys.getenv("prot")
    gene <- Sys.getenv("gene")
    chr <- Sys.getenv("chr")
    region37 <- Sys.getenv("region37")
    HPC_WORK <- Sys.getenv("HPC_WORK")
    f <- file.path(find.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
    chain <- rtracklayer::import.chain(f)
    gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
    vcf <- file.path("~/INF/METAL/gwas2vcf",paste0(prot,".vcf.gz"))
    gwas_stats <- gwasvcf::query_gwas(vcf, chrompos = region37) %>%
                  gwasvcf::vcf_to_granges() %>%
                  keepSeqlevels(chr) %>%
                  renameSeqlevels(paste0("chr",chr))
    gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
      unlist() %>%
      dplyr::as_tibble() %>%
      dplyr::transmute(chromosome = seqnames,
                       position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
      dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
      dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(row_count = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(row_count == 1) %>%
      mutate(chromosome=gsub("chr","",chromosome))
    save(gwas_stats_hg38,file=paste0("~/INF/coloc/sumstats/",gene,".rda"))
  '
}

function eQTL()
{
  export src=~/rds/public_databases/GTEx/csv
  if [ ! -d ${INF}/coloc/GTEx ]; then mkdir -p ${INF}/coloc/GTEx; fi
  ls ${src}/*gz | \
  parallel -C' ' '
    export tissue=$(echo {} | xargs -l basename -s .tsv.gz)
    cat <(head -1 ~/pQTLtools/inst/eQTL-Catalogue/column_names.GTEx) \
        <(tabix {} ${region} | grep ${ensGene}) > ${INF}/coloc/GTEx/${gene}-${tissue}.tsv
  '
  export src=~/rds/public_databases/eQTLCatalogue
  if [ ! -d ${INF}/coloc/ge ]; then mkdir -p ${INF}/coloc/ge; fi
  ls ${src}/*gz | \
  parallel -C' ' '
    export tissue=$(echo {} | xargs -l basename -s .tsv.gz)
    cat <(head -1 ~/pQTLtools/inst/eQTL-Catalogue/column_names.Alasoo) \
        <(tabix {} ${region} | grep ${ensGene})> ${INF}/coloc/ge/${gene}-${tissue}.tsv
  '
}

function run_coloc()
{
  if [ ! -d ${INF}/coloc/coloc ]; then mkdir -p ${INF}/coloc/coloc; fi
  for src in GTEx ge
  do
    export src=${src}
    ls ${INF}/coloc/${src}/${gene}* | \
    grep -e CCL23 -e CCL25 -e CCL4 -e CD6 -e FGF5 -e CCL3 -e TNFSF14 | \
    parallel -C' ' '
      export tissue=$(echo {} | xargs -l basename -s .tsv)
      Rscript -e "
        suppressMessages(library(dplyr))
        gene <- Sys.getenv(\"gene\")
        load(paste0(\"~/INF/coloc/sumstats/\",gene,\".rda\"))
        gwas_stats <- mutate(gwas_stats_hg38,snpid=gap::chr_pos_a1_a2(chromosome,position,REF,ALT))
        src <- Sys.getenv(\"src\")
        tissue <- Sys.getenv(\"tissue\")
        eqtl_stats <- read.delim(paste0(\"~/INF/coloc/\",src,\"/\",tissue,\".tsv\")) %>%
                      mutate(snpid=gap::chr_pos_a1_a2(chromosome,position,ref,alt)) %>%
                      left_join(select(gwas_stats,snpid,REF,ALT),by=\"snpid\") %>%
                      mutate(maf=as.numeric(maf),sign=if_else(alt==ALT,1,-1),beta=sign*beta) %>%
                      filter(!is.na(beta))
        dup1 <- duplicated(with(gwas_stats,snpid))
        dup2 <- duplicated(with(eqtl_stats,snpid))
        run_coloc <- function(gwas_sumstats=gwas_stats[!dup1,],eqtl_sumstats=eqtl_stats[!dup2,])
        {
          eQTL_dataset <- with(eqtl_sumstats, list(beta=beta,varbeta=se^2,N=an,MAF=maf,type=\"quant\",snp=snpid))
          gwas_dataset <- with(gwas_sumstats, list(beta=ES,varbeta=SE^2,type=\"quant\",snp=snpid,MAF=MAF,N=SS))
          coloc_res <- coloc::coloc.abf(dataset1=gwas_dataset, dataset2=eQTL_dataset, p1=1e-4, p2=1e-4, p12=1e-5)
          res_formatted <- dplyr::as_tibble(t(as.data.frame(coloc_res$summary)))
        }
        res <- run_coloc()
        r <- Sys.getenv(\"SLURM_ARRAY_TASK_ID\")
        write.table(res,file=file.path(\"~/INF/coloc/coloc\",paste0(r,\"-\",src,\"-\",tissue,\".out\")))
      "
    '
  done
}

function run()
{
  cis
  sumstats38
  eQTL
  run_coloc
}

run_coloc

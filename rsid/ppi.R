# Identification of potential drug targets
# Using PPI networks for alternative drug targets search

library("magrittr")
library("dplyr")
library("purrr")
library("glue")
library("epigraphdb")

GENE_NAME <- "IL12B"
OUTCOME_TRAIT <- "Multiple sclerosis"

get_drug_targets_ppi <- function(gene_name) {
  endpoint <- "/gene/druggability/ppi"
  params <- list(gene_name = gene_name)
  df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  df
}

ppi_df <- get_drug_targets_ppi(gene_name = GENE_NAME)
ppi_df

# gene list

get_gene_list <- function(ppi_df, include_primary_gene = TRUE) {
  if (include_primary_gene) {
    gene_list <- c(
      ppi_df %>% pull(`g1.name`) %>% unique(),
      ppi_df %>% filter(`g2.druggability_tier` == "Tier 1") %>% pull(`g2.name`)
    )
  } else {
    gene_list <- ppi_df %>%
      filter(`g2.druggability_tier` == "Tier 1") %>%
      pull(`g2.name`)
  }
  gene_list
}

gene_list <- get_gene_list(ppi_df)
gene_list

# Using Mendelian randomization results for causal effect estimation

extract_mr <- function(outcome_trait, gene_list, qtl_type) {
  endpoint <- "/xqtl/single-snp-mr"
  per_gene <- function(gene_name) {
    params <- list(
      exposure_gene = gene_name,
      outcome_trait = outcome_trait,
      qtl_type = qtl_type,
      pval_threshold = 1e-5
    )
    df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
    df
  }
  res_df <- gene_list %>% map_df(per_gene)
  res_df
}

xqtl_df <- c("pQTL", "eQTL") %>% map_df(function(qtl_type) {
  extract_mr(
    outcome_trait = OUTCOME_TRAIT,
    gene_list = gene_list,
    qtl_type = qtl_type
  ) %>%
    mutate(qtl_type = qtl_type)
})
xqtl_df

# Using literature evidence for results enrichment and triangulation

extract_literature <- function(outcome_trait, gene_list) {
  per_gene <- function(gene_name) {
    endpoint <- "/gene/literature"
    params <- list(
      gene_name = gene_name,
      object_name = outcome_trait %>% stringr::str_to_lower()
    )
    df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
    df
  }
  res_df <- gene_list %>% map_df(per_gene)
  res_df %>%
    mutate(literature_count = map_int(pubmed_id, function(x) length(x)))
}

literature_df <- extract_literature(
  outcome_trait = OUTCOME_TRAIT,
  gene_list = gene_list
)
literature_df

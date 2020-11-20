# Distinguishing vertical and horizontal pleiotropy for SNP-protein associations

options(width=200)

# SCALLOP/INF data

m <- read.delim("work/INF1.merge")
head(m)
dim(m)
mp <-with(m,table(MarkerName,prot))
apply(mp,2,sum)
length(apply(mp,2,sum))
subset(m,prot=="SCF")
subset(m,prot=="uPA")
length(apply(mp,1,sum))
mp <- apply(mp,1,sum)
mp[mp>4]
subset(m,MarkerName%in%c("chr1:159175354_A_G", "chr12:111884608_C_T"))
subset(m,MarkerName%in%c("chr1:159175354_A_G"))
subset(m,MarkerName%in%("chr12:111884608_C_T"))
library(gap)
subset(inf1,prot%in%subset(m,MarkerName%in%("chr1:159175354_A_G"))$prot)$gene
subset(inf1,prot%in%subset(m,MarkerName%in%("chr12:111884608_C_T"))$prot)$gene

library("magrittr")
library("dplyr")
library("purrr")
library("glue")
library("igraph")
library("epigraphdb")

PPI_N_INTERMEDIATE_PROTEINS <- 1
# chr1:159175354_A_G
SNP <- "rs12075"
GENELIST <- c("CXCL11","CXCL10","CD5","IL12B","CXCL9","CD244")

# chr12:111884608_C_T
SNP <- "rs3184504"
GENELIST <- c("CCL2","CCL11","CCL8","CCL7","CXCL6","CCL13")

# Mapping genes to proteins
get_gene_protein <- function(genelist) {
  endpoint <- "/mappings/gene-to-protein"
  params <- list(
    gene_name_list = genelist %>% I()
  )
  r <- query_epigraphdb(
    route = endpoint,
    params = params,
    mode = "table",
    method = "POST"
  )
  protein_df <- r
  if (nrow(protein_df) > 0) {
    res_df <- protein_df %>%
      select(protein_name = `gene.name`, uniprot_id = `protein.uniprot_id`)
  } else {
    res_df <- tibble() %>% set_names(c("protein_name", "uniprot_id"))
  }
  res_df
}

gene_protein_df <- get_gene_protein(genelist = GENELIST)
gene_protein_df

# Identifying the involvement in the same biological pathways
# For each protein we retrieve the pathways they are involved in
get_protein_pathway <- function(gene_protein_df) {
  endpoint <- "/protein/in-pathway"
  params <- list(
    uniprot_id_list = gene_protein_df %>% pull(`uniprot_id`) %>% I()
  )
  df <- query_epigraphdb(route = endpoint, params = params, mode = "table", method = "POST")

  if (nrow(df) > 0) {
    res_df <- gene_protein_df %>%
      select(`uniprot_id`) %>%
      left_join(df, by = c("uniprot_id"))
  } else {
    res_df <- gene_protein_df %>%
      select(`uniprot_id`) %>%
      mutate(pathway_count = NA_integer_, pathway_reactome_id = NA_character_)
  }
  res_df <- res_df %>%
    mutate(
      pathway_count = ifelse(is.na(pathway_count), 0L, as.integer(pathway_count)),
      pathway_reactome_id = ifelse(is.na(pathway_reactome_id), c(), pathway_reactome_id)
    )
  res_df
}
pathway_df <- get_protein_pathway(gene_protein_df = gene_protein_df)
pathway_df

# for each pair of proteins we match the pathways they have in common
get_shared_pathway <- function(pathway_df) {
  # For the protein-pathway data
  # Get protein-protein permutations where they share pathways

  per_permutation <- function(pathway_df, permutation) {
    df <- pathway_df %>% filter(uniprot_id %in% permutation)
    primary_pathway <- pathway_df %>%
      filter(uniprot_id == permutation[1]) %>%
      pull(pathway_reactome_id) %>%
      unlist()
    assoc_pathway <- pathway_df %>%
      filter(uniprot_id == permutation[2]) %>%
      pull(pathway_reactome_id) %>%
      unlist()
    intersect(primary_pathway, assoc_pathway)
  }

  pairwise_permutations <- pathway_df %>%
    pull(`uniprot_id`) %>%
    gtools::permutations(n = length(.), r = 2, v = .)
  shared_pathway_df <- tibble(
    protein = pairwise_permutations[, 1],
    assoc_protein = pairwise_permutations[, 2]
  ) %>%
    mutate(
      shared_pathway = map2(`protein`, `assoc_protein`, function(x, y) {
        per_permutation(pathway_df = pathway_df, permutation = c(x, y))
      }),
      combination = map2_chr(`protein`, `assoc_protein`, function(x, y) {
        comb <- sort(c(x, y))
        paste(comb, collapse = ",")
      }),
      count = map_int(`shared_pathway`, function(x) na.omit(x) %>% length()),
      connected = count > 0
    )
  shared_pathway_df
}
shared_pathway_df <- get_shared_pathway(pathway_df)
n_pairs <- length(shared_pathway_df %>% filter(count > 0))
print(glue::glue("Num. shared_pathway pairs: {n_pairs}"))
#> Num. shared_pathway pairs: 6
shared_pathway_df %>% arrange(desc(count))

# detailed pathway information
get_pathway_info <- function(reactome_id) {
  endpoint <- "/meta/nodes/Pathway/search"
  params <- list(id = reactome_id)
  df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  df
}

pathway <- shared_pathway_df %>%
  pull(shared_pathway) %>%
  unlist() %>%
  unique()

pathway_info <- pathway %>% map_df(get_pathway_info)
pathway_info %>% print()

# count the number of nodes in each connected community and plot the graph.

protein_df_to_graph <- function(df) {
  df_connected <- df %>%
    filter(connected) %>%
    distinct(`combination`, .keep_all = TRUE)
  nodes <- df %>%
    pull(protein) %>%
    unique()
  graph <- igraph::graph_from_data_frame(
    df_connected,
    directed = FALSE, vertices = nodes
  )
  graph$layout <- igraph::layout_with_kk
  graph
}

graph_to_protein_groups <- function(graph) {
  graph %>%
    igraph::components() %>%
    igraph::groups() %>%
    tibble(group_member = .) %>%
    mutate(group_size = map_int(`group_member`, length)) %>%
    arrange(desc(group_size))
}

pathway_protein_graph <- shared_pathway_df %>% protein_df_to_graph()
pathway_protein_groups <- pathway_protein_graph %>% graph_to_protein_groups()
pathway_protein_groups %>% str()

# Retrieving shared protein-protein interactions

get_ppi <- function(gene_protein_df, n_intermediate_proteins = 0) {
  endpoint <- "/protein/ppi/pairwise"
  params <- list(
    uniprot_id_list = gene_protein_df %>% pull("uniprot_id") %>% I(),
    n_intermediate_proteins = n_intermediate_proteins
  )
  df <- query_epigraphdb(
    route = endpoint,
    params = params,
    mode = "table",
    method = "POST"
  )

  if (nrow(df) > 0) {
    res_df <- gene_protein_df %>%
      select(protein = uniprot_id) %>%
      left_join(df, by = c("protein"))
  } else {
    res_df <- gene_protein_df %>%
      select(protein = uniprot_id) %>%
      mutate(assoc_protein = NA_character_, path_size = NA_integer_)
  }
  res_df
}
ppi_df <- get_ppi(gene_protein_df, n_intermediate_proteins = PPI_N_INTERMEDIATE_PROTEINS)
ppi_df

get_shared_ppi <- function(ppi_df) {
  per_row <- function(ppi_df, x, y) {
    res <- ppi_df %>%
      filter(
        `protein` == x,
        `assoc_protein` == y
      ) %>%
      nrow() %>%
      `>`(0)
    res
  }

  pairwise_permutations <- ppi_df %>%
    pull(`protein`) %>%
    unique() %>%
    gtools::permutations(n = length(.), r = 2, v = .)
  shared_ppi_df <- tibble(
    protein = pairwise_permutations[, 1],
    assoc_protein = pairwise_permutations[, 2]
  ) %>%
    mutate(
      shared_ppi = map2_lgl(`protein`, `assoc_protein`, function(x, y) {
        per_row(ppi_df = ppi_df, x = x, y = y)
      }),
      combination = map2_chr(`protein`, `assoc_protein`, function(x, y) {
        comb <- sort(c(x, y))
        paste(comb, collapse = ",")
      }),
      connected = shared_ppi
    )
  shared_ppi_df
}

shared_ppi_df <- get_shared_ppi(ppi_df)
n_ppi_pairs <- shared_ppi_df %>%
  filter(shared_ppi) %>%
  nrow()
print(glue::glue("Num. shared_ppi pairs: {n_ppi_pairs}"))
#> Num. shared_ppi pairs: 10
shared_ppi_df

ppi_protein_graph <- protein_df_to_graph(shared_ppi_df)
ppi_protein_groups <- graph_to_protein_groups(ppi_protein_graph)
ppi_protein_groups %>% str()
plot(ppi_protein_graph)

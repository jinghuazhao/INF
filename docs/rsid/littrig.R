## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("magrittr")
library("dplyr")
library("purrr")
library("glue")
library("ggplot2")
library("igraph")
library("epigraphdb")

## -----------------------------------------------------------------------------
STARTING_TRAIT <- "Inflammatory bowel disease"

## -----------------------------------------------------------------------------
get_mr <- function(trait) {
  endpoint <- "/mr"
  params <- list(
    exposure_trait = trait,
    pval_threshold = 1e-10
  )
  mr_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  mr_df
}

mr_df <- get_mr(STARTING_TRAIT)
mr_df %>% glimpse()

## -----------------------------------------------------------------------------
trait_to_disease <- function(trait) {
  endpoint <- "/ontology/gwas-efo-disease"
  params <- list(trait = trait)
  disease_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  if (nrow(disease_df) > 0) {
    res <- disease_df %>% pull(`disease.label`)
  } else {
    res <- c()
  }
  res
}

disease_df <- mr_df %>%
  mutate(disease = map(`outcome.trait`, trait_to_disease)) %>%
  filter(map_lgl(`disease`, function(x) !is.null(x)))
disease_df

## -----------------------------------------------------------------------------
get_gwas_pair_literature <- function(gwas_id, assoc_gwas_id) {
  endpoint <- "/literature/gwas/pairwise"
  # NOTE in this example we blacklist to semmentic types
  params <- list(
    gwas_id = gwas_id,
    assoc_gwas_id = assoc_gwas_id,
    by_gwas_id = TRUE,
    pval_threshold = 1e-1,
    semmantic_types = "nusq",
    semmantic_types = "dsyn",
    blacklist = TRUE,
    limit = 1000
  )
  lit_df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  lit_df
}

GWAS_ID_X <- "ieu-a-1088"
GWAS_ID_Y <- "ieu-a-6"
lit_df <- get_gwas_pair_literature(GWAS_ID_X, GWAS_ID_Y)

glimpse(lit_df)

## -----------------------------------------------------------------------------
# Predicate counts for SemMed triples for trait X
lit_df %>%
  count(`s1.predicate`) %>%
  arrange(desc(n))

## -----------------------------------------------------------------------------
# Predicate counts for SemMed triples for trait Y
lit_df %>%
  count(`s2.predicate`) %>%
  arrange(desc(n))

## -----------------------------------------------------------------------------
# Filter out some predicates that are not informative
pred_filter <- c("COEXISTS_WITH", "ASSOCIATED_WITH")
lit_df_filter <- lit_df %>%
  filter(
    !`s1.predicate` %in% pred_filter,
    !`s2.predicate` %in% pred_filter
  )
lit_df_filter %>%
  count(`s1.predicate`) %>%
  arrange(desc(n))

## -----------------------------------------------------------------------------
lit_df_filter %>%
  count(`s2.predicate`) %>%
  arrange(desc(n))

## -----------------------------------------------------------------------------
lit_counts <- lit_df_filter %>%
  count(`st.type`, `st.name`) %>%
  arrange(`st.type`, desc(`n`))
lit_counts %>% print(n = 30)

## -----------------------------------------------------------------------------
lit_counts %>%
  filter(n < 100) %>%
  {
    ggplot(.) +
      aes(x = `st.name`, y = `n`, fill = `st.type`) +
      geom_col() +
      geom_text(
        aes(label = `n`),
        position = position_dodge(0.9),
        hjust = 0
      ) +
      coord_flip()
  }

## -----------------------------------------------------------------------------
focus_term <- "Leptin"
lit_detail <- lit_df_filter %>% filter(`st.name` == focus_term)
lit_detail %>% head()

## -----------------------------------------------------------------------------
lit_detail <- lit_detail %>%
  mutate_at(vars(`gwas.trait`, `assoc_gwas.trait`), stringr::str_to_upper)

# add node types: 1 - selected GWAS, 2 - traits from literature, 3 - current focus term connecting 1 and 2
nodes <- bind_rows(
  lit_detail %>% select(node = `gwas.trait`) %>% distinct() %>% mutate(node_type = 1),
  lit_detail %>% select(node = `assoc_gwas.trait`) %>% distinct() %>% mutate(node_type = 1),
  lit_detail %>% select(node = `s1.subject_name`) %>% distinct() %>% mutate(node_type = 2),
  lit_detail %>% select(node = `s2.subject_name`) %>% distinct() %>% mutate(node_type = 2),
  lit_detail %>% select(node = `s1.object_name`) %>% distinct() %>% mutate(node_type = 2),
  lit_detail %>% select(node = `s2.object_name`) %>% distinct() %>% mutate(node_type = 2)
) %>%
  mutate(node_type = ifelse(node == focus_term, 3, node_type)) %>%
  distinct()
nodes

## -----------------------------------------------------------------------------
edges <- bind_rows(
  # exposure -> s1 subject
  lit_detail %>%
    select(node = `gwas.trait`, assoc_node = `s1.subject_name`) %>%
    distinct(),
  # s2 object -> outcome
  lit_detail %>%
    select(node = `s2.object_name`, assoc_node = `assoc_gwas.trait`) %>%
    distinct(),
  # s1 subject - s1 predicate -> s1 object
  lit_detail %>%
    select(
      node = `s1.subject_name`, assoc_node = `s1.object_name`,
      label = `s1.predicate`
    ) %>%
    distinct(),
  # s2 subject - s2 predicate -> s2 object
  lit_detail %>%
    select(
      node = `s2.subject_name`, assoc_node = `s2.object_name`,
      label = `s2.predicate`
    ) %>%
    distinct()
) %>%
  distinct()
edges

## -----------------------------------------------------------------------------
plot_network <- function(edges, nodes, show_edge_labels = FALSE) {

  # default is to not display edge labels
  if (!show_edge_labels) {
    edges <- select(edges, -label)
  }

  graph <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  graph$layout <- layout_with_kk

  # generate colors based on node type
  colrs <- c("tomato", "steelblue", "gold")
  V(graph)$color <- colrs[V(graph)$node_type]

  # Configure canvas
  default_mar <- par("mar")
  new_mar <- c(0, 0, 0, 0)
  par(mar = new_mar)

  plot.igraph(graph,
    vertex.size = 13,
    vertex.label.color = "black",
    vertex.label.family = "Helvetica",
    vertex.label.cex = 0.8,
    edge.arrow.size = 0.4,
    edge.label.color = "black",
    edge.label.family = "Helvetica",
    edge.label.cex = 0.5
  )
  par(mar = default_mar)
}

## ----dpi=300------------------------------------------------------------------
plot_network(edges, nodes, show_edge_labels = TRUE)

## -----------------------------------------------------------------------------
get_literature <- function(gwas_id, semmed_triple_id) {
  endpoint <- "/literature/gwas"
  params <- list(
    gwas_id = gwas_id,
    semmed_triple_id = semmed_triple_id,
    by_gwas_id = TRUE,
    pval_threshold = 1e-1
  )
  df <- query_epigraphdb(route = endpoint, params = params, mode = "table")
  df %>% select(`triple.id`, `lit.pubmed_id`)
}

pub_df <- bind_rows(
  lit_detail %>%
    select(gwas_id = `gwas.id`, semmed_triple_id = `s1.id`) %>%
    distinct(),
  lit_detail %>%
    select(gwas_id = `assoc_gwas.id`, semmed_triple_id = `s2.id`) %>%
    distinct()
) %>%
  transpose() %>%
  map_df(function(x) get_literature(x$gwas_id, x$semmed_triple_id))
pub_df

## -----------------------------------------------------------------------------
sessionInfo()


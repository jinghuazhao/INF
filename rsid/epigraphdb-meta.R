## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("magrittr")
library("dplyr")
library("purrr")
library("igraph")
library("epigraphdb")

## -----------------------------------------------------------------------------
endpoint <- "/meta/schema"
params <- list(
  graphviz = FALSE,
  plot = FALSE
)
metadata <- query_epigraphdb(
  route = endpoint, params = params, mode = "raw"
)

metadata %>% str(2)

## -----------------------------------------------------------------------------
meta_node_df <- metadata %>%
  pluck("nodes") %>%
  {
    names <- names(.)
    transpose(.) %>%
      as_tibble() %>%
      mutate(meta_node = names) %>%
      # Hide properties column which does not display well
      select(meta_node, count) %>%
      # We also need to flatten count
      mutate(count = flatten_int(count))
  }

meta_node_df %>%
  arrange(meta_node) %>%
  mutate(count = format(count, big.mark = ","))

## -----------------------------------------------------------------------------
meta_rel_df <- metadata %>%
  pluck("edges") %>%
  {
    names <- names(.)
    transpose(.) %>%
      as_tibble() %>%
      mutate(meta_rel = names) %>%
      mutate(count = flatten_int(count)) %>%
      select(meta_rel, count)
  } %>%
  inner_join(
    metadata %>% pluck("connections") %>% {
      transpose(.) %>%
        as_tibble() %>%
        mutate(meta_rel = flatten_chr(rel)) %>%
        mutate_at(vars(from_node, to_node), flatten_chr) %>%
        select(meta_rel, from_node, to_node)
    }
  )

meta_rel_df %>%
  arrange(from_node, to_node) %>%
  mutate(count = format(count, big.mark = ","))

## -----------------------------------------------------------------------------
graph <- meta_rel_df %>%
  select(from_node, to_node) %>%
  igraph::graph_from_data_frame(directed = FALSE)
graph$layout <- igraph::layout_with_kk
plot(graph)

## -----------------------------------------------------------------------------
endpoint <- "/meta/nodes/id-name-schema"
meta_node_fields <- query_epigraphdb(
  route = endpoint, params = NULL, mode = "raw"
)
meta_node_fields

## -----------------------------------------------------------------------------
name <- "body mass index"

endpoint <- "/meta/nodes/Gwas/search"
params <- list(name = name)

results <- query_epigraphdb(
  route = endpoint, params = params, mode = "table"
)
results

## -----------------------------------------------------------------------------
id <- "ieu-a-2"

endpoint <- "/meta/nodes/Gwas/search"
params <- list(id = id)

results <- query_epigraphdb(
  route = endpoint, params = params, mode = "table"
)
results

## -----------------------------------------------------------------------------
query <- "
    MATCH (exposure:Gwas)-[mr:MR]->(outcome:Gwas)
    WHERE exposure.trait = 'Body mass index'
    RETURN exposure, outcome, mr LIMIT 2
"

endpoint <- "/cypher"
params <- list(query = query)

# NOTE this is a POST request
results <- query_epigraphdb(
  route = endpoint, params = params, method = "POST",
  mode = "table"
)
results

## -----------------------------------------------------------------------------
sessionInfo()


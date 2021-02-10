# ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("epigraphdb")
url <- getOption("epigraphdb.api.url")
payload <- list(
  exposure_trait = "Body mass index",
  outcome_trait = "Coronary heart disease"
)
r <- httr::RETRY("GET", glue::glue("{url}/mr"), query = payload)
r %>%
  httr::content(as = "parsed") %>%
  str(2)

## ----epigraphdb-mr------------------------------------------------------------
epigraphdb::mr(
  exposure_trait = "Body mass index",
  outcome_trait = "Coronary heart disease",
  mode = "raw"
) %>%
  str(2)

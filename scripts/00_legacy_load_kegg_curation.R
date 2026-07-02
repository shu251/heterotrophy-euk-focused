library(tidyverse)

kegg_curation <- read_csv("data/raw/kegg_curation_SKHU_May2026.csv")

# Coverage of the heterotrophy classification scheme
kegg_curation %>%
  count(CLASSIFICATION_BROAD, sort = TRUE)

kegg_curation %>%
  filter(is.na(CLASSIFICATION_BROAD) | CLASSIFICATION_BROAD == "") %>%
  count(CLASS, sort = TRUE)

# Rows still missing a PFAM cross-reference
kegg_curation %>%
  summarise(
    n_total = n(),
    n_missing_pfam = sum(is.na(PFAM) | PFAM == "")
  )

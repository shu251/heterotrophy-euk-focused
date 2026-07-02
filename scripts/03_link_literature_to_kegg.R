library(tidyverse)
library(fuzzyjoin)

lit <- read_csv("data/raw/literature_heterotrophy.csv", show_col_types = FALSE)
kegg <- read_csv("data/processed/kegg_classified_all.csv", show_col_types = FALSE)

# --- Direct join for rows that already have ko_id ---
direct <- lit |>
  filter(!is.na(ko_id)) |>
  left_join(
    kegg |> select(KEGG, NAME_SHORT, NAME_FULL, CLASSIFICATION_BROAD, CLASSIFICATION_SPECIFIC, DESCRIPTION),
    by = c("ko_id" = "KEGG")
  ) |>
  mutate(ko_match_confidence = "exact")

# --- Fuzzy join for rows without ko_id ---
needs_match <- lit |> filter(is.na(ko_id))

kegg_lookup <- kegg |>
  mutate(search_text = paste(coalesce(NAME_SHORT, ""), coalesce(NAME_FULL, ""), sep = " ") |> tolower()) |>
  select(KEGG, NAME_SHORT, NAME_FULL, CLASSIFICATION_BROAD, CLASSIFICATION_SPECIFIC, DESCRIPTION, search_text)

# Fuzzy join on gene_name against KEGG name fields (Jaro-Winkler >= 0.9)
if (nrow(needs_match) > 0 && requireNamespace("stringdist", quietly = TRUE)) {
  fuzzy_matched <- stringdist_left_join(
    needs_match |> mutate(query = tolower(gene_name)),
    kegg_lookup,
    by = c("query" = "search_text"),
    method = "jw",
    max_dist = 0.1,
    distance_col = "jw_dist"
  ) |>
    group_by(ref_id, gene_name) |>
    slice_min(jw_dist, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      ko_id = coalesce(ko_id, KEGG),
      ko_match_confidence = if_else(!is.na(KEGG), "fuzzy", "unmatched")
    ) |>
    select(-query, -search_text, -jw_dist)
} else {
  fuzzy_matched <- needs_match |>
    mutate(
      KEGG = NA_character_, NAME_SHORT = NA_character_, NAME_FULL = NA_character_,
      CLASSIFICATION_BROAD = NA_character_, CLASSIFICATION_SPECIFIC = NA_character_,
      DESCRIPTION = NA_character_, ko_match_confidence = "unmatched"
    )
}

result <- bind_rows(direct, fuzzy_matched) |>
  rename(
    kegg_name_short = NAME_SHORT,
    kegg_name_full = NAME_FULL,
    kegg_classification_broad = CLASSIFICATION_BROAD,
    kegg_classification_specific = CLASSIFICATION_SPECIFIC,
    kegg_description = DESCRIPTION
  ) |>
  arrange(function_broad, gene_name)

write_csv(result, "data/processed/lit_kegg_linked.csv")

cat("=== KEGG Linkage Summary ===\n")
result |> count(ko_match_confidence) |> print()
cat("\nSaved: data/processed/lit_kegg_linked.csv (", nrow(result), "rows)\n")

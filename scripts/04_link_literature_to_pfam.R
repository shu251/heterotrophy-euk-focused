library(tidyverse)
library(fuzzyjoin)

lit <- read_csv("data/raw/literature_heterotrophy.csv", show_col_types = FALSE)
pfam <- read_csv("data/raw/pfam_parsed_table.csv", show_col_types = FALSE)

# Standardize PFAM column names (inspect what's in the file)
cat("PFAM table columns:", paste(names(pfam), collapse = ", "), "\n")

# Direct join for rows with pfam_id provided
direct <- lit |>
  filter(!is.na(pfam_id)) |>
  left_join(pfam, by = c("pfam_id" = names(pfam)[1])) |>
  mutate(pfam_match_confidence = "exact")

# Fuzzy join gene_name against PFAM description for unlinked rows
needs_match <- lit |> filter(is.na(pfam_id))

# Find the description column
desc_col <- names(pfam)[str_detect(tolower(names(pfam)), "desc|name")]
if (length(desc_col) == 0) desc_col <- names(pfam)[2]

pfam_lookup <- pfam |>
  mutate(search_text = tolower(.data[[desc_col[1]]])) |>
  select(1, all_of(desc_col[1]), search_text)

if (nrow(needs_match) > 0 && requireNamespace("stringdist", quietly = TRUE)) {
  fuzzy_matched <- stringdist_left_join(
    needs_match |> mutate(query = tolower(gene_name)),
    pfam_lookup,
    by = c("query" = "search_text"),
    method = "jw",
    max_dist = 0.15,
    distance_col = "jw_dist"
  ) |>
    group_by(ref_id, gene_name) |>
    slice_min(jw_dist, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      pfam_id = coalesce(pfam_id, as.character(.data[[names(pfam_lookup)[1]]])),
      pfam_match_confidence = if_else(!is.na(.data[[names(pfam_lookup)[1]]]), "fuzzy", "unmatched")
    ) |>
    select(-query, -search_text, -jw_dist)
} else {
  fuzzy_matched <- needs_match |>
    mutate(pfam_match_confidence = "unmatched")
}

result <- bind_rows(direct, fuzzy_matched) |>
  arrange(function_broad, gene_name)

write_csv(result, "data/processed/lit_pfam_linked.csv")

cat("=== PFAM Linkage Summary ===\n")
result |> count(pfam_match_confidence) |> print()
cat("\nSaved: data/processed/lit_pfam_linked.csv (", nrow(result), "rows)\n")

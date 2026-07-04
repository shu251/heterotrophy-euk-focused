library(tidyverse)
library(fuzzyjoin)

lit <- read_csv("data/raw/literature_heterotrophy.csv", show_col_types = FALSE)
pfam <- read_csv("data/raw/pfam_parsed_table.csv", show_col_types = FALSE)

# Standardize PFAM column names (inspect what's in the file)
cat("PFAM table columns:", paste(names(pfam), collapse = ", "), "\n")

# Identify the PFAM ID column (Pfam_ID) and description column
pfam_id_col  <- names(pfam)[str_detect(names(pfam), regex("pfam_id", ignore_case = TRUE))][1]
if (is.na(pfam_id_col)) pfam_id_col <- "Pfam_ID"
# Prefer the plain "Description" column for fuzzy matching (not "Descriptive_ID")
desc_col <- if ("Description" %in% names(pfam)) "Description" else {
  candidates <- names(pfam)[str_detect(tolower(names(pfam)), "^desc")]
  if (length(candidates) > 0) candidates[1] else names(pfam)[2]
}
cat("Using PFAM ID col:", pfam_id_col, " | Description col:", desc_col, "\n")

# Direct join for rows with pfam_id provided
direct <- lit |>
  filter(!is.na(pfam_id)) |>
  left_join(pfam, by = c("pfam_id" = pfam_id_col)) |>
  mutate(pfam_match_confidence = "exact")

# Fuzzy join gene_name against PFAM description for unlinked rows
needs_match <- lit |> filter(is.na(pfam_id))

pfam_lookup <- pfam |>
  mutate(search_text = tolower(.data[[desc_col]])) |>
  select(all_of(pfam_id_col), all_of(desc_col), search_text)

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
      pfam_id = coalesce(pfam_id, as.character(.data[[pfam_id_col]])),
      pfam_match_confidence = if_else(!is.na(.data[[pfam_id_col]]), "fuzzy", "unmatched")
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

library(tidyverse)

lit <- read_csv("data/raw/literature_heterotrophy.csv", show_col_types = FALSE)

# --- Validation checks ---
required_cols <- c("ref_id", "gene_name", "function_broad")
missing_required <- lit |>
  filter(if_any(all_of(required_cols), is.na)) |>
  select(all_of(required_cols))

if (nrow(missing_required) > 0) {
  message("WARNING: ", nrow(missing_required), " rows missing required fields:")
  print(missing_required)
} else {
  message("OK: All rows have required fields (ref_id, gene_name, function_broad)")
}

dupes <- lit |>
  count(ref_id, gene_name) |>
  filter(n > 1)

if (nrow(dupes) > 0) {
  message("WARNING: ", nrow(dupes), " duplicate ref_id + gene_name combinations:")
  print(dupes)
} else {
  message("OK: No duplicate ref_id + gene_name combinations")
}

# --- Summary report ---
cat("\n=== Literature Input Summary ===\n")
cat("Total rows:", nrow(lit), "\n")
cat("Unique references:", n_distinct(lit$ref_id), "\n")
cat("Unique genes:", n_distinct(lit$gene_name), "\n")
cat("Rows with ko_id:", sum(!is.na(lit$ko_id)), "\n")
cat("Rows with pfam_id:", sum(!is.na(lit$pfam_id)), "\n")
cat("Rows with cazy_family:", sum(!is.na(lit$cazy_family)), "\n\n")

cat("Breakdown by function_broad:\n")
lit |> count(function_broad, sort = TRUE) |> print()

# --- Save summary ---
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

gene_ref_summary <- lit |>
  group_by(gene_name, gene_alias, function_broad, function_specific) |>
  summarise(
    n_refs = n(),
    ref_ids = paste(ref_id, collapse = "; "),
    ko_id = first(na.omit(ko_id)),
    pfam_id = first(na.omit(pfam_id)),
    cazy_family = first(na.omit(cazy_family)),
    .groups = "drop"
  ) |>
  arrange(function_broad, gene_name)

write_csv(gene_ref_summary, "data/processed/gene_ref_summary.csv")
message("\nSaved: data/processed/gene_ref_summary.csv (", nrow(gene_ref_summary), " unique genes)")

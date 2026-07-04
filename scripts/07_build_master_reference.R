library(tidyverse)

# Load all processed outputs
lit_kegg <- read_csv("data/processed/lit_kegg_linked.csv", show_col_types = FALSE)
lit_pfam_raw <- read_csv("data/processed/lit_pfam_linked.csv", show_col_types = FALSE)

# Build a clean pfam lookup: one row per ref_id × gene_name with the best pfam_id
# lit_pfam has the original pfam_id col + a Pfam_ID col from the join — coalesce them
lit_pfam <- lit_pfam_raw |>
  mutate(
    pfam_id_resolved = coalesce(
      if ("Pfam_ID" %in% names(lit_pfam_raw)) .data[["Pfam_ID"]] else NA_character_,
      pfam_id
    )
  ) |>
  select(ref_id, gene_name, pfam_id_resolved, pfam_match_confidence) |>
  rename(pfam_id_pfam = pfam_id_resolved, pfam_match_conf = pfam_match_confidence)

# CAZy reference: rows that have gene_name linking back to literature
cazy_raw <- read_csv("data/processed/cazy_reference.csv", show_col_types = FALSE)
cazy <- cazy_raw |>
  filter(heterotrophy_relevant) |>
  mutate(gene_name = as.character(gene_name), ref_id = as.character(ref_id)) |>
  filter(!is.na(gene_name)) |>
  select(family_id, class, description, gene_name, function_broad, ref_id) |>
  rename(cazy_family_match = family_id, cazy_class = class, cazy_description = description)

# Add PFAM linkage — resolve pfam_id: prefer literature-supplied, then PFAM-matched
master <- lit_kegg |>
  left_join(lit_pfam, by = c("ref_id", "gene_name")) |>
  mutate(
    pfam_id = coalesce(pfam_id, pfam_id_pfam)
  ) |>
  select(-pfam_id_pfam)

# Add CAZy family info where literature entries specify one
master <- master |>
  left_join(
    cazy |> select(gene_name, ref_id, cazy_family_match, cazy_class, cazy_description),
    by = c("gene_name", "ref_id"),
    relationship = "many-to-many"
  ) |>
  mutate(cazy_family = coalesce(cazy_family, cazy_family_match)) |>
  select(-cazy_family_match)

# Collapse to one row per gene × function, refs as semicolon list
collapsed <- master |>
  group_by(gene_name, gene_alias, function_broad, function_specific, ko_id,
           kegg_name_short, kegg_name_full, kegg_classification_broad,
           pfam_id, cazy_family,
           cazy_class = if ("cazy_class" %in% names(master)) cazy_class else NA_character_) |>
  summarise(
    ref_ids = paste(sort(unique(ref_id)), collapse = "; "),
    dois = paste(sort(unique(doi)), collapse = "; "),
    years = paste(sort(unique(year)), collapse = "; "),
    organism_studied = paste(sort(unique(na.omit(organism_studied))), collapse = "; "),
    ko_match_confidence = first(na.omit(ko_match_confidence)),
    n_refs = n_distinct(ref_id),
    notes = paste(unique(na.omit(notes)), collapse = " | "),
    .groups = "drop"
  ) |>
  arrange(function_broad, function_specific, gene_name)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(collapsed, "data/processed/heterotrophy_master_reference.csv")

cat("=== Master Reference Summary ===\n")
cat("Total unique gene × function entries:", nrow(collapsed), "\n")
cat("Entries with KEGG KO:", sum(!is.na(collapsed$ko_id)), "\n")
cat("Entries with PFAM:", sum(!is.na(collapsed$pfam_id)), "\n")
cat("Entries with CAZy:", sum(!is.na(collapsed$cazy_family)), "\n\n")

cat("Breakdown by function_broad:\n")
collapsed |> count(function_broad, sort = TRUE) |> print()

cat("\nTop entries by number of supporting references:\n")
collapsed |>
  arrange(desc(n_refs)) |>
  select(gene_name, function_broad, ko_id, n_refs, ref_ids) |>
  head(10) |>
  print()

message("\nSaved: data/processed/heterotrophy_master_reference.csv")

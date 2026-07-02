library(tidyverse)

kegg_raw <- read_csv("data/raw/kegg_curation_SKHU_May2026.csv", show_col_types = FALSE)

# Keyword patterns by heterotrophy category
patterns <- list(
  Phagocytosis = c("phago", "lysosom", "vacuol", "engulf", "endocyt", "phagolyso",
                   "phagosome", "endosome", "macropinocyt"),
  `Vesicle trafficking` = c("clathrin", "dynamin", "\\brab\\b", "\\bvps\\b", "escrt",
                             "snare", "\\bvamp\\b", "\\bsec\\d", "retromer", "ap-1",
                             "ap-2", "ap-3", "coat protein", "coatomer", "arf\\b",
                             "copi", "copii", "sorting nexin"),
  Digestion = c("lipase", "protease", "peptidase", "amylase", "glucosidase",
                "cathepsin", "acid hydrolase", "lysosomal", "hydrolase",
                "chitinase", "nuclease", "phospholipase", "esterase",
                "beta-galactosidase", "beta-glucuronidase", "alpha-mannosidase"),
  Assimilation = c("beta.oxidation", "fatty.acid.*degrad", "acetyl.coa",
                   "slc.*transport", "amino acid transport", "sugar transport",
                   "pufa", "peroxisom"),
  `Prey detection` = c("lectin", "chemoreceptor", "flagell", "ciliar", "motil"),
  `Motility/cytoskeletal` = c("\\bactin\\b", "myosin", "kinesin", "dynein",
                               "tubulin", "cytoskelet", "filopod", "lamellopod")
)

search_text <- function(df) {
  paste(
    coalesce(df$NAME_FULL, ""),
    coalesce(df$NAME_SHORT, ""),
    coalesce(df$DESCRIPTION, ""),
    sep = " "
  ) |> tolower()
}

classified <- kegg_raw |>
  filter(!is.na(CLASSIFICATION_BROAD) & CLASSIFICATION_BROAD != "")

unclassified <- kegg_raw |>
  filter(is.na(CLASSIFICATION_BROAD) | CLASSIFICATION_BROAD == "")

text <- search_text(unclassified)

newly_classified <- unclassified |> mutate(text = text)

for (cat_name in names(patterns)) {
  pat <- paste(patterns[[cat_name]], collapse = "|")
  hit <- str_detect(newly_classified$text, pat)
  newly_classified <- newly_classified |>
    mutate(
      CLASSIFICATION_BROAD = case_when(
        hit & is.na(CLASSIFICATION_BROAD) ~ cat_name,
        hit & CLASSIFICATION_BROAD == "" ~ cat_name,
        TRUE ~ CLASSIFICATION_BROAD
      ),
      CLASSIFICATION_SPECIFIC = case_when(
        hit & (is.na(CLASSIFICATION_SPECIFIC) | CLASSIFICATION_SPECIFIC == "") ~ cat_name,
        TRUE ~ CLASSIFICATION_SPECIFIC
      )
    )
  rm(hit)
}

newly_classified <- newly_classified |>
  filter(!is.na(CLASSIFICATION_BROAD) & CLASSIFICATION_BROAD != "") |>
  mutate(classification_source = "keyword_auto") |>
  select(-text)

combined <- bind_rows(
  classified |> mutate(classification_source = "manual"),
  newly_classified
)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(combined, "data/processed/kegg_classified_all.csv")

cat("=== KEGG Classification Summary ===\n")
cat("Originally classified:", nrow(classified), "\n")
cat("Newly classified from unclassified pool:", nrow(newly_classified), "\n")
cat("Total classified:", nrow(combined), "\n")
cat("Still unclassified:", nrow(unclassified) - nrow(newly_classified), "\n\n")
cat("Breakdown by CLASSIFICATION_BROAD:\n")
combined |> count(CLASSIFICATION_BROAD, sort = TRUE) |> print()

message("\nSaved: data/processed/kegg_classified_all.csv")

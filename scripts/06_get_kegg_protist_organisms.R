library(tidyverse)
library(httr2)

# KEGG REST API: list all organisms
message("Downloading KEGG organism list...")
resp <- request("https://rest.kegg.jp/list/organism") |>
  req_timeout(60) |>
  req_perform()

raw_text <- resp_body_string(resp)

orgs <- read_tsv(
  I(raw_text),
  col_names = c("T_number", "org_code", "name", "taxonomy"),
  show_col_types = FALSE
)

cat("Total KEGG organisms:", nrow(orgs), "\n")

# Protist lineage filter terms
protist_terms <- c(
  # Stramenopiles
  "Stramenopiles", "Ochrophyta", "Oomycota", "Labyrinthulomycetes",
  "Diatom", "Bacillariophyta", "Chrysophyceae", "Phaeophyceae", "Xanthophyceae",
  # Alveolata
  "Apicomplexa", "Plasmodium", "Toxoplasma", "Cryptosporidium", "Theileria", "Babesia",
  "Ciliophora", "Tetrahymena", "Paramecium", "Stentor",
  "Dinoflagellata", "Symbiodinium", "Symbiodiniaceae",
  # Rhizaria
  "Rhizaria", "Foraminifera", "Radiolaria", "Chlorarachniophyta",
  # Excavata
  "Euglenozoa", "Trypanosoma", "Leishmania", "Euglena",
  "Diplomonadida", "Giardia",
  "Parabasalia", "Trichomonas",
  "Heterolobosea", "Naegleria",
  # Amoebozoa
  "Amoebozoa", "Dictyostelium", "Entamoeba", "Acanthamoeba",
  # Choanoflagellatea
  "Choanoflagellatea", "Monosiga", "Salpingoeca",
  # Haptophyta
  "Haptophyta", "Emiliania", "Prymnesium",
  # Cryptophyta
  "Cryptophyta", "Cryptomonas", "Guillardia",
  # Chlorophyta (unicellular / basal)
  "Ostreococcus", "Micromonas", "Bathycoccus", "Chlamy"
)

pattern <- paste(protist_terms, collapse = "|")

protists <- orgs |>
  filter(str_detect(taxonomy, regex(pattern, ignore_case = TRUE)) |
         str_detect(name, regex(pattern, ignore_case = TRUE))) |>
  arrange(taxonomy, name)

cat("Protist organisms identified:", nrow(protists), "\n")
cat("\nBreakdown by major group:\n")

protists |>
  mutate(major_group = case_when(
    str_detect(taxonomy, regex("Stramenopiles|Ochrophyta|Oomycota|Labyrinthulomycetes|Diatom|Bacillariophyta", TRUE)) ~ "Stramenopiles",
    str_detect(taxonomy, regex("Apicomplexa|Plasmodium|Toxoplasma|Cryptosporidium|Theileria|Babesia", TRUE)) ~ "Apicomplexa",
    str_detect(taxonomy, regex("Ciliophora|Tetrahymena|Paramecium", TRUE)) ~ "Ciliophora",
    str_detect(taxonomy, regex("Dinoflagellata|Symbiodinium", TRUE)) ~ "Dinoflagellata",
    str_detect(taxonomy, regex("Euglenozoa|Trypanosoma|Leishmania|Euglena", TRUE)) ~ "Euglenozoa",
    str_detect(taxonomy, regex("Diplomonadida|Giardia", TRUE)) ~ "Diplomonadida",
    str_detect(taxonomy, regex("Amoebozoa|Dictyostelium|Entamoeba", TRUE)) ~ "Amoebozoa",
    str_detect(taxonomy, regex("Choanoflagellatea|Monosiga|Salpingoeca", TRUE)) ~ "Choanoflagellatea",
    str_detect(taxonomy, regex("Haptophyta|Emiliania|Prymnesium", TRUE)) ~ "Haptophyta",
    str_detect(taxonomy, regex("Cryptophyta|Cryptomonas", TRUE)) ~ "Cryptophyta",
    TRUE ~ "Other/Unclassified protist"
  )) |>
  count(major_group, sort = TRUE) |>
  print()

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(protists, "data/processed/kegg_protist_organisms.csv")
message("\nSaved: data/processed/kegg_protist_organisms.csv")

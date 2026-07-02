library(tidyverse)
library(httr2)
library(rvest)

# CAZy class pages to scrape
cazy_classes <- list(
  GH  = list(url = "https://www.cazy.org/Glycoside-Hydrolases.html",  heterotrophy_relevant = TRUE,  priority = "high"),
  PL  = list(url = "https://www.cazy.org/Polysaccharide-Lyases.html", heterotrophy_relevant = TRUE,  priority = "high"),
  CE  = list(url = "https://www.cazy.org/Carbohydrate-Esterases.html",heterotrophy_relevant = TRUE,  priority = "high"),
  AA  = list(url = "https://www.cazy.org/Auxiliary-Activities.html",   heterotrophy_relevant = TRUE,  priority = "medium"),
  GT  = list(url = "https://www.cazy.org/GlycosylTransferases.html",   heterotrophy_relevant = FALSE, priority = "low"),
  CBM = list(url = "https://www.cazy.org/Carbohydrate-Binding-Modules.html", heterotrophy_relevant = FALSE, priority = "low")
)

parse_cazy_page <- function(class_abbrev, class_info) {
  message("Fetching CAZy class: ", class_abbrev)
  Sys.sleep(1)  # polite delay

  tryCatch({
    page <- read_html(class_info$url)

    # Family links follow pattern like "GH1", "GH2", etc.
    links <- page |> html_elements("a") |>
      html_attr("href") |>
      na.omit()

    family_links <- links[str_detect(links, paste0("/", class_abbrev, "\\d+\\.html$"))]
    family_ids <- str_extract(family_links, paste0(class_abbrev, "\\d+"))

    # Get table rows with family + description
    tables <- page |> html_elements("table")
    rows <- tibble(family_id = character(), description = character())

    for (tbl in tables) {
      tbl_data <- tryCatch(html_table(tbl, fill = TRUE), error = function(e) NULL)
      if (!is.null(tbl_data) && ncol(tbl_data) >= 2) {
        col1 <- tbl_data[[1]]
        col2 <- tbl_data[[2]]
        hits <- str_detect(as.character(col1), paste0("^", class_abbrev, "\\d+$"))
        if (any(hits, na.rm = TRUE)) {
          rows <- bind_rows(rows, tibble(
            family_id = as.character(col1[hits]),
            description = as.character(col2[hits])
          ))
        }
      }
    }

    if (nrow(rows) == 0 && length(family_ids) > 0) {
      rows <- tibble(family_id = unique(family_ids), description = NA_character_)
    }

    rows |> mutate(
      class = class_abbrev,
      heterotrophy_relevant = class_info$heterotrophy_relevant,
      priority = class_info$priority
    )
  }, error = function(e) {
    message("  Failed to fetch ", class_abbrev, ": ", conditionMessage(e))
    tibble(family_id = NA_character_, description = NA_character_,
           class = class_abbrev,
           heterotrophy_relevant = class_info$heterotrophy_relevant,
           priority = class_info$priority)
  })
}

cazy_all <- map2_dfr(names(cazy_classes), cazy_classes, parse_cazy_page) |>
  filter(!is.na(family_id)) |>
  distinct()

# Join to any CAZy families already cited in literature
lit <- read_csv("data/raw/literature_heterotrophy.csv", show_col_types = FALSE)
lit_cazy <- lit |>
  filter(!is.na(cazy_family)) |>
  select(cazy_family, gene_name, function_broad, ref_id) |>
  rename(family_id = cazy_family)

cazy_final <- cazy_all |>
  left_join(lit_cazy, by = "family_id", relationship = "many-to-many") |>
  arrange(class, family_id)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(cazy_final, "data/processed/cazy_reference.csv")

cat("=== CAZy Family Summary ===\n")
cazy_all |> count(class, heterotrophy_relevant, priority) |> print()
cat("\nTotal families downloaded:", nrow(cazy_all), "\n")
cat("Saved: data/processed/cazy_reference.csv\n")

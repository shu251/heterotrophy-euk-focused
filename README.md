# Heterotrophy Transcripts — Euk-Focused

Reference annotation framework for identifying heterotrophy-related functional
traits (phagotrophy, osmotrophy, parasitism, saprotrophy) in protist
metatranscriptomes, by cross-referencing **KEGG**, **PFAM**, and **CAZy**
enzyme/family databases.

## Data (`data/raw/`)

Copied from [`../KEGG_DB`](../KEGG_DB) on 2026-07-01. That folder remains the
working source for KEGG/PFAM curation across projects — treat files here as a
snapshot for this project, not the canonical copy.

| File | Rows | Description |
|---|---|---|
| `kegg_curation_SKHU_May2026.csv` | 1,545 | Primary curated KO list. Columns: `CLASSIFICATION_BROAD`, `CLASSIFICATION_SPECIFIC`, `CLASS`, `KEGG`, `GO`, `PFAM`, `NAME_SHORT`, `NAME_FULL`, `DESCRIPTION`, `REFERENCE`, `WHY`. Most complete/recent curation. |
| `trophy_curated_kegg.csv` | 455 | Earlier, simpler KO list (KEGG ID + `CLASS_REFINE` only, e.g. "Fatty acid breakdown"). Likely superseded by the May2026 file — kept for reference/cross-check. |
| `reformat-kegg-pfam-skh.csv` | 475 | Earlier KEGG+PFAM curation pass (carbon fixation, gluconeogenesis, etc.). |
| `pfam_parsed_table.csv` | 19,632 | Full PFAM-A family reference table (Pfam ID, description, gathering threshold, clan). Not yet linked to the heterotrophy KO list. |

No CAZy data exists yet in any prior project — that leg needs to be built from
scratch (e.g. via the [dbCAN](https://bcb.unl.edu/dbCAN2/) family list or the
[CAZy website](http://www.cazy.org/)).

## Status (2026-07-01)

- **KEGG**: 1,545 KOs curated. Classified so far: Phagotrophy (300),
  Parasitism (72), Nutrient processing (47), Cellular function (21).
  **1,105 rows still have blank `CLASSIFICATION_BROAD`/`CLASSIFICATION_SPECIFIC`**
  (mostly cytoskeletal/trafficking genes pulled from Labarre et al. 2019 —
  e.g. dynein, clathrin — not yet placed into the heterotrophy scheme).
- **PFAM**: reference table present, not cross-linked to curated KOs.
- **CAZy**: not started.

## Next steps

1. Finish classifying the ~1,105 blank rows in `kegg_curation_SKHU_May2026.csv`.
2. Link `PFAM` column in the KO curation against `pfam_parsed_table.csv`.
3. Build a CAZy family list and cross-reference against carbohydrate-active
   heterotrophy traits.
4. Identify KEGG organism IDs for target protists
   (https://rest.kegg.jp/list/organism) to scope which KOs are even present
   in the taxa of interest.

## Structure

- `data/raw/` — source curation files (see table above)
- `data/processed/` — cleaned/merged outputs
- `scripts/` — R scripts / Quarto docs for loading, merging, and analyzing data
- `output/` — figures, summary tables

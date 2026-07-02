# Heterotrophy Transcripts — Euk-Focused

Reference annotation framework for identifying heterotrophy-related genes, proteins,
and pathways in protist metatranscriptomes. **Literature-first**: every entry traces
back to a published reference. KEGG KOs, PFAM domains, and CAZy enzyme families
serve as the cross-reference layer.

## Architecture

```
literature_heterotrophy.csv  ← primary input (human-editable; add new refs here)
         ↓
01_validate  →  02_classify_kegg  →  03_link_kegg  →  04_link_pfam
                                                            ↓
                                         05_get_cazy  →  06_get_kegg_orgs
                                                            ↓
                                              07_build_master_reference.csv
```

## Data (`data/raw/`)

| File | Description |
|------|-------------|
| `literature_heterotrophy.csv` | **Primary input** — one row per gene per reference. Edit this to add new papers. See column guide below. |
| `kegg_curation_SKHU_May2026.csv` | 1,545 curated KOs (snapshot from `../KEGG_DB`); 440 classified, 1,105 pending. |
| `pfam_parsed_table.csv` | 19,632 PFAM-A families (reference table). |
| `trophy_curated_kegg.csv` | Earlier KO list (455 rows); kept for cross-check. |
| `reformat-kegg-pfam-skh.csv` | Earlier KEGG+PFAM pass; kept for cross-check. |

### `literature_heterotrophy.csv` columns

| Column | Required | Description |
|--------|----------|-------------|
| `ref_id` | ✓ | Short unique key, e.g. `Labarre2020` |
| `doi` | | DOI string |
| `first_author` | | Last name |
| `year` | | Publication year |
| `journal` | | Journal name |
| `title` | | Paper title |
| `gene_name` | ✓ | Gene/protein name as used in the paper |
| `gene_alias` | | Semicolon-separated aliases |
| `function_broad` | ✓ | e.g. `Phagocytosis`, `Digestion`, `Vesicle trafficking`, `Assimilation`, `Detection`, `Scaling up` |
| `function_specific` | | e.g. `Phagosome formation`, `Lysosomal acidification` |
| `organism_studied` | | Protist taxon(a) studied |
| `ko_id` | | KEGG KO ID (K#####); script 03 fills blanks via fuzzy match |
| `pfam_id` | | PFAM accession (PF#####); script 04 fills blanks |
| `cazy_family` | | CAZy family (e.g. `GH18`); script 05 links |
| `notes` | | Free text |

**To add a new reference:** append rows (one per gene per paper) and re-run scripts 01 → 07.

## Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `01_validate_literature_input.R` | Validate CSV; print summary by function_broad |
| `02_classify_kegg_unclassified.R` | Keyword-classify 1,105 blank KEGG KOs |
| `03_link_literature_to_kegg.R` | Join lit genes → KEGG KOs (exact + fuzzy) |
| `04_link_literature_to_pfam.R` | Join lit genes → PFAM families |
| `05_get_cazy_families.R` | Download CAZy GH/PL/CE/AA families from cazy.org |
| `06_get_kegg_protist_organisms.R` | Fetch KEGG organism codes for protists |
| `07_build_master_reference.R` | Merge all → `heterotrophy_master_reference.csv` |

Run in order: `01` → `02` → `03` → `04` → `05` → `06` → `07`

## Outputs (`data/processed/`)

| File | Description |
|------|-------------|
| `gene_ref_summary.csv` | One row per unique gene; ref count, KEGG/PFAM/CAZy IDs |
| `kegg_classified_all.csv` | All KEGG KOs with classification (manual + keyword-auto) |
| `lit_kegg_linked.csv` | Literature table + KEGG KO matches |
| `lit_pfam_linked.csv` | Literature table + PFAM matches |
| `cazy_reference.csv` | CAZy GH/PL/CE/AA family table |
| `kegg_protist_organisms.csv` | KEGG organism codes for protists |
| `heterotrophy_master_reference.csv` | **Final output** — one row per gene × function |

## Key KEGG pathway anchors
| KEGG ID | Pathway | Heterotrophy role |
|---------|---------|-------------------|
| ko04145 | Phagosome | Prey engulfment |
| ko04142 | Lysosome | Intracellular digestion |
| ko00071 | Fatty acid degradation | Prey lipid assimilation |

## Seed literature (10 core references, 2026-07-02)
Pernthaler (2005) · Caron et al. (2017) · Kiørboe (2024) · Labarre et al. (2020, 2021) ·
Massana et al. (2021) · Wootton et al. (2007) · Boulais et al. (2010) ·
Lie et al. (2022) · Carradec et al. (2018)

## MCP Server (interactive use)

A local MCP server lets you add references and search the database from any Claude conversation.

**Setup** (one time):
```bash
cd mcp
# venv already created at mcp/.venv — install is done
# Register the server by pointing Claude Code at .mcp.json in the project root
```

**Available tools:**
| Tool | Description |
|------|-------------|
| `add_reference(doi, gene_name, function_broad, ...)` | Auto-fetches paper metadata from CrossRef; appends validated row to CSV |
| `search_genes(query)` | Fuzzy search gene_name + gene_alias |
| `get_coverage(function_broad?)` | Summary counts by function category |
| `link_gene_to_kegg(gene_name)` | Live KEGG REST lookup for candidate KO IDs |

**Environment variable:** `HETEROTROPHY_CSV` — override default CSV path (set in `.mcp.json`).

## Test Dataset

`data/raw/test_literature_heterotrophy.csv` — 5 rows, 2 refs, 5 genes with known KO/PFAM/CAZy IDs.
Use to validate scripts 01–07 before running the full 100-row dataset.
Expected: Digestion=3, Phagocytosis=1, Vesicle trafficking=1; all 4 genes with ko_id exact-matched in script 03.

## Structure
- `data/raw/` — source files (KEGG, PFAM, literature CSV, test CSV)
- `data/processed/` — pipeline outputs
- `scripts/` — R scripts (run 01–07 in order)
- `mcp/` — Python MCP server for interactive reference addition
- `output/` — figures, summary tables

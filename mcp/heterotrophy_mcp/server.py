"""
Heterotrophy Reference MCP Server

Tools:
  add_reference    — guided row entry with CrossRef DOI auto-fill
  search_genes     — fuzzy search across gene_name + gene_alias
  get_coverage     — summary counts by function_broad
  link_gene_to_kegg — live KEGG REST lookup for candidate KO IDs
"""

import csv
import json
import os
import re
from difflib import SequenceMatcher
from pathlib import Path

import httpx
from mcp.server.fastmcp import FastMCP

CSV_PATH = Path(
    os.environ.get(
        "HETEROTROPHY_CSV",
        Path(__file__).parent.parent.parent / "data" / "raw" / "literature_heterotrophy.csv",
    )
)

COLUMNS = [
    "ref_id", "doi", "first_author", "year", "journal", "title",
    "gene_name", "gene_alias", "function_broad", "function_specific",
    "organism_studied", "ko_id", "pfam_id", "cazy_family", "notes",
]

FUNCTION_BROAD_OPTIONS = [
    "Phagocytosis", "Vesicle trafficking", "Digestion", "Assimilation",
    "Detection", "Osmotrophy", "Egestion", "Scaling up",
]

mcp = FastMCP("heterotrophy-reference")


def _read_rows() -> list[dict]:
    if not CSV_PATH.exists():
        return []
    with open(CSV_PATH, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _append_row(row: dict) -> None:
    exists = CSV_PATH.exists()
    with open(CSV_PATH, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        if not exists:
            writer.writeheader()
        writer.writerow({col: row.get(col, "") for col in COLUMNS})


def _fetch_crossref(doi: str) -> dict:
    url = f"https://api.crossref.org/works/{doi.strip()}"
    try:
        r = httpx.get(url, timeout=10, headers={"User-Agent": "heterotrophy-mcp/1.0 (mailto:skhu@tamu.edu)"})
        r.raise_for_status()
        item = r.json()["message"]
        authors = item.get("author", [])
        first_author = authors[0].get("family", "") if authors else ""
        year = ""
        if "published" in item:
            parts = item["published"].get("date-parts", [[]])
            year = str(parts[0][0]) if parts and parts[0] else ""
        elif "published-print" in item:
            parts = item["published-print"].get("date-parts", [[]])
            year = str(parts[0][0]) if parts and parts[0] else ""
        journal = ""
        if "container-title" in item and item["container-title"]:
            journal = item["container-title"][0]
        title = ""
        if "title" in item and item["title"]:
            title = item["title"][0]
        return {"first_author": first_author, "year": year, "journal": journal, "title": title}
    except Exception as e:
        return {"error": str(e)}


def _similarity(a: str, b: str) -> float:
    return SequenceMatcher(None, a.lower(), b.lower()).ratio()


@mcp.tool()
def add_reference(
    doi: str,
    gene_name: str,
    function_broad: str,
    gene_alias: str = "",
    function_specific: str = "",
    organism_studied: str = "",
    ko_id: str = "",
    pfam_id: str = "",
    cazy_family: str = "",
    notes: str = "",
) -> str:
    """
    Add a gene entry to the heterotrophy reference CSV.
    Fetches paper metadata (author, year, journal, title) from CrossRef using the DOI.
    Required: doi, gene_name, function_broad.
    function_broad must be one of: Phagocytosis, Vesicle trafficking, Digestion,
    Assimilation, Detection, Osmotrophy, Egestion, Scaling up.
    """
    if not doi or not gene_name or not function_broad:
        return "Error: doi, gene_name, and function_broad are all required."

    if function_broad not in FUNCTION_BROAD_OPTIONS:
        return (
            f"Error: function_broad '{function_broad}' not recognized. "
            f"Choose one of: {', '.join(FUNCTION_BROAD_OPTIONS)}"
        )

    meta = _fetch_crossref(doi)
    if "error" in meta:
        return f"CrossRef lookup failed ({meta['error']}). Check the DOI and try again."

    # Build ref_id as AuthorYear
    ref_id = re.sub(r"[^A-Za-z0-9]", "", meta["first_author"]) + meta["year"]
    if not ref_id:
        ref_id = doi.split("/")[-1][:12]

    # Check for duplicate
    existing = _read_rows()
    for row in existing:
        if row.get("ref_id") == ref_id and row.get("gene_name", "").lower() == gene_name.lower():
            return f"Entry already exists: {ref_id} / {gene_name}. No row added."

    row = {
        "ref_id": ref_id,
        "doi": doi,
        "first_author": meta["first_author"],
        "year": meta["year"],
        "journal": meta["journal"],
        "title": meta["title"],
        "gene_name": gene_name,
        "gene_alias": gene_alias,
        "function_broad": function_broad,
        "function_specific": function_specific,
        "organism_studied": organism_studied,
        "ko_id": ko_id,
        "pfam_id": pfam_id,
        "cazy_family": cazy_family,
        "notes": notes,
    }
    _append_row(row)
    return (
        f"Added: {ref_id} / {gene_name} ({function_broad})\n"
        f"  Paper: {meta['first_author']} {meta['year']} — {meta['journal']}\n"
        f"  Title: {meta['title']}\n"
        f"  KO: {ko_id or '—'}  PFAM: {pfam_id or '—'}  CAZy: {cazy_family or '—'}"
    )


@mcp.tool()
def search_genes(query: str, top_n: int = 10) -> str:
    """
    Search the literature CSV for genes matching a query string.
    Searches gene_name and gene_alias. Returns top matches with ref_ids and function.
    """
    rows = _read_rows()
    if not rows:
        return "CSV is empty or not found."

    scored = []
    q = query.lower()
    for row in rows:
        name = row.get("gene_name", "").lower()
        alias = row.get("gene_alias", "").lower()
        score = max(
            _similarity(q, name),
            max((_similarity(q, a.strip()) for a in alias.split(";") if a.strip()), default=0),
            1.0 if q in name or q in alias else 0,
        )
        scored.append((score, row))

    scored.sort(key=lambda x: -x[0])
    seen = set()
    results = []
    for score, row in scored[:top_n * 3]:
        if score < 0.3:
            break
        key = (row["gene_name"], row["ref_id"])
        if key in seen:
            continue
        seen.add(key)
        results.append(
            f"  {row['gene_name']} ({row.get('gene_alias','')}) | {row['function_broad']} | "
            f"ref: {row['ref_id']} | KO: {row.get('ko_id','—')} | PFAM: {row.get('pfam_id','—')}"
        )
        if len(results) >= top_n:
            break

    if not results:
        return f"No matches found for '{query}'."
    return f"Results for '{query}':\n" + "\n".join(results)


@mcp.tool()
def get_coverage(function_broad: str = "") -> str:
    """
    Summarize how many genes and references are in the literature CSV,
    optionally filtered to one function_broad category.
    """
    rows = _read_rows()
    if not rows:
        return "CSV is empty or not found."

    if function_broad:
        rows = [r for r in rows if r.get("function_broad", "").lower() == function_broad.lower()]
        if not rows:
            return f"No entries found for function_broad='{function_broad}'."

    from collections import Counter
    by_func = Counter(r.get("function_broad", "") for r in rows)
    n_genes = len(set(r.get("gene_name", "") for r in rows))
    n_refs = len(set(r.get("ref_id", "") for r in rows))
    n_ko = sum(1 for r in rows if r.get("ko_id"))
    n_pfam = sum(1 for r in rows if r.get("pfam_id"))
    n_cazy = sum(1 for r in rows if r.get("cazy_family"))

    lines = [
        f"Total rows: {len(rows)}",
        f"Unique genes: {n_genes}  |  Unique refs: {n_refs}",
        f"Has KO ID: {n_ko}  |  Has PFAM: {n_pfam}  |  Has CAZy: {n_cazy}",
        "",
        "By function_broad:",
    ]
    for func, count in by_func.most_common():
        lines.append(f"  {func}: {count}")
    return "\n".join(lines)


@mcp.tool()
def link_gene_to_kegg(gene_name: str) -> str:
    """
    Query the KEGG REST API for a gene name and return candidate KO IDs.
    Useful for looking up the KO for a gene before adding it to the CSV.
    """
    url = f"https://rest.kegg.jp/find/ko/{gene_name}"
    try:
        r = httpx.get(url, timeout=10)
        r.raise_for_status()
        lines = [l for l in r.text.strip().splitlines() if l]
        if not lines:
            return f"No KEGG KO matches found for '{gene_name}'."
        results = []
        for line in lines[:10]:
            parts = line.split("\t", 1)
            ko = parts[0].replace("ko:", "")
            desc = parts[1] if len(parts) > 1 else ""
            results.append(f"  {ko}  {desc}")
        return f"KEGG KO candidates for '{gene_name}':\n" + "\n".join(results)
    except Exception as e:
        return f"KEGG query failed: {e}"


if __name__ == "__main__":
    mcp.run()

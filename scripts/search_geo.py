#!/usr/bin/env python3
"""Search NCBI GEO for omics datasets via Entrez E-utilities."""

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET

from scripts.utils import fetch_with_retry, resolve_doi, setup_logger
from scripts.sample_utils import build_sample_table, SAMPLE_CAP
from scripts.supplement_fetch import fetch_supplementary_tables

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

CHARACTERISTIC_KEY_MAP = {
    "tissue_site": ["tissue", "organ", "body site", "anatomic site"],
    "sample_type": ["tumor type", "sample type", "primary", "metastasis", "metastatic"],
    "treatment_status": ["treatment", "therapy", "chemo", "drug"],
    "disease": ["disease", "diagnosis", "pathology", "histology", "disease state"],
    "cell_type": ["cell type", "cell line", "sorted population"],
    "age": ["age"],
    "sex": ["sex", "gender"],
    "stage": ["stage", "grade", "tnm"],
}

OMIC_TO_GEO_TERMS = {
    "scRNAseq": '("single cell" OR "scRNA" OR "10x" OR "Smart-seq") AND "Expression profiling by high throughput sequencing"[Filter]',
    "bulkRNAseq": '"Expression profiling by high throughput sequencing"[Filter] NOT "single cell" NOT "scRNA"',
    "scMultiomicseq": '("single cell" OR "multiome") AND ("ATAC" OR "multi-omic")',
    "scGenomicseq": '("single cell" OR "single-cell") AND ("DNA sequencing" OR "whole genome")',
    "bulkGenomicseq": '("whole genome sequencing" OR "WGS" OR "exome") NOT "single cell"',
    "ATACseq": '("ATAC-seq" OR "ATACseq" OR "chromatin accessibility")',
    "ChIPseq": '("ChIP-seq" OR "ChIPseq" OR "chromatin immunoprecipitation sequencing")',
    "CRISPR": '("CRISPR" OR "genome-wide screen" OR "loss of function screen")',
}

logger = setup_logger("search_geo")


def build_esearch_query(omic, organism="human", disease=None, tissue=None, date_from=None, date_to=None):
    parts = []
    omic_term = OMIC_TO_GEO_TERMS.get(omic, f'"{omic}"')
    parts.append(omic_term)
    organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
    org = organism_map.get(organism.lower(), organism)
    parts.append(f'"{org}"[Organism]')
    if disease:
        parts.append(f'"{disease}"')
    if tissue:
        parts.append(f'"{tissue}"')
    if date_from or date_to:
        df = date_from.replace("-", "/") if date_from else "1900/01/01"
        dt = date_to.replace("-", "/") if date_to else "3000/12/31"
        parts.append(f'"{df}"[PDAT] : "{dt}"[PDAT]')
    return " AND ".join(parts)


def _get_item_text(doc, name):
    for item in doc.findall("Item"):
        if item.get("Name") == name:
            return (item.text or "").strip()
    return ""


def _get_item_int(doc, name):
    text = _get_item_text(doc, name)
    try:
        return int(text)
    except (ValueError, TypeError):
        return 0


def _get_pubmed_ids(doc):
    for item in doc.findall("Item"):
        if item.get("Name") == "PubMedIds":
            return [sub.text.strip() for sub in item.findall("Item") if sub.text]
    return []


def _infer_data_files(supp_file):
    files = []
    if not supp_file:
        return [{"type": "supplementary", "format": "unknown", "size_mb": 0}]
    for fmt in supp_file.split(";"):
        fmt = fmt.strip().upper()
        if fmt:
            files.append({"type": "count_matrix" if fmt in ("H5", "MTX", "CSV", "TSV") else "supplementary", "format": fmt.lower(), "size_mb": 0})
    return files if files else [{"type": "supplementary", "format": "unknown", "size_mb": 0}]


def _assess_metadata_quality(doc):
    has_title = bool(_get_item_text(doc, "title"))
    has_summary = bool(_get_item_text(doc, "summary"))
    has_samples = _get_item_int(doc, "n_samples") > 0
    has_pubmed = bool(_get_pubmed_ids(doc))
    score = sum([has_title, has_summary, has_samples, has_pubmed])
    if score >= 4:
        return "good"
    elif score >= 2:
        return "partial"
    return "minimal"


def parse_geo_record(doc, omic):
    accession = _get_item_text(doc, "Accession")
    if not accession:
        accession = f"GSE{_get_item_text(doc, 'GSE')}" if _get_item_text(doc, "GSE") else ""
    pdat = _get_item_text(doc, "PDAT").replace("/", "-")
    summary = _get_item_text(doc, "summary")
    pubmed_ids = _get_pubmed_ids(doc)
    paper = {"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": summary, "pubmed_id": pubmed_ids[0] if pubmed_ids else ""}
    return {
        "accession": accession, "source": "GEO", "title": _get_item_text(doc, "title"),
        "organism": _get_item_text(doc, "taxon"), "omic_type": omic,
        "platform": _get_item_text(doc, "GPL") or "unknown",
        "disease": "", "tissue": "", "sample_count": _get_item_int(doc, "n_samples"),
        "condition_groups": [], "data_files": _infer_data_files(_get_item_text(doc, "suppFile")),
        "paper": paper, "metadata_quality": _assess_metadata_quality(doc), "date_submitted": pdat,
    }


def _map_characteristic(key: str, value: str):
    """Map a SOFT characteristic key/value to a (schema_column, value) tuple, or None.

    Matching strategy: first try exact equality against each keyword, then fall
    back to substring containment.  The two-pass approach prevents shorter
    keywords (e.g. "age") from shadowing longer ones (e.g. "stage") when the
    short keyword happens to be a substring of the key.
    """
    key_lower = key.lower().strip()

    # Pass 1: exact match
    for schema_col, keywords in CHARACTERISTIC_KEY_MAP.items():
        for kw in keywords:
            if key_lower == kw:
                return (schema_col, value.strip())

    # Pass 2: substring containment — process entries with longer keywords first
    # so that e.g. "stage" (5 chars) is checked before "age" (3 chars).
    for schema_col, keywords in sorted(CHARACTERISTIC_KEY_MAP.items(), key=lambda item: max(len(kw) for kw in item[1]), reverse=True):
        for kw in keywords:
            if kw in key_lower:
                return (schema_col, value.strip())

    return None


def parse_soft_samples(soft_text: str) -> list[dict]:
    """Parse GEO SOFT format text into sample metadata rows.

    Splits on ``^SAMPLE`` blocks, extracts ``!Sample_characteristics_ch1``
    key:value pairs (mapped via CHARACTERISTIC_KEY_MAP), and
    ``!Sample_platform_id``. Caps at SAMPLE_CAP.
    """
    rows = []
    current: dict | None = None

    for line in soft_text.splitlines():
        line = line.rstrip()
        if line.startswith("^SAMPLE"):
            if current is not None:
                rows.append(current)
                if len(rows) >= SAMPLE_CAP:
                    break
            # Start new sample block
            parts = line.split("=", 1)
            sample_id = parts[1].strip() if len(parts) == 2 else ""
            current = {"sample_id": sample_id}
        elif current is not None:
            if line.startswith("!Sample_characteristics_ch1"):
                # Format: !Sample_characteristics_ch1 = key: value
                rest = line.split("=", 1)[1].strip() if "=" in line else ""
                if ":" in rest:
                    char_key, char_val = rest.split(":", 1)
                    mapped = _map_characteristic(char_key, char_val)
                    if mapped is not None:
                        schema_col, val = mapped
                        if schema_col not in current:
                            current[schema_col] = val
            elif line.startswith("!Sample_platform_id"):
                rest = line.split("=", 1)[1].strip() if "=" in line else ""
                current["platform_id"] = rest

    # Append last block if within cap
    if current is not None and len(rows) < SAMPLE_CAP:
        rows.append(current)

    return rows


def fetch_geo_samples(accession: str, paper: dict | None = None) -> dict | None:
    """Fetch SOFT for *accession*, parse samples, optionally enrich from supplementary tables.

    Returns a build_sample_table dict (source="GEO") or None on failure.
    """
    try:
        params = {"acc": accession, "form": "text", "view": "brief"}
        resp = fetch_with_retry(GEO_SOFT_URL, params=params)
        rows = parse_soft_samples(resp.text)
    except Exception as exc:
        logger.warning("fetch_geo_samples failed for %s: %s", accession, exc)
        return None

    # Attempt supplementary enrichment when paper metadata is available
    if paper and (paper.get("pubmed_id") or paper.get("doi")):
        sample_ids = [r["sample_id"] for r in rows if r.get("sample_id")]
        try:
            enriched = fetch_supplementary_tables(
                doi=paper.get("doi") or None,
                pmid=paper.get("pubmed_id") or None,
                sample_ids=sample_ids or None,
            )
            if enriched:
                id_to_row = {r["sample_id"]: r for r in rows}
                for enr in enriched:
                    sid = enr.get("sample_id")
                    if sid and sid in id_to_row:
                        for k, v in enr.items():
                            if k != "sample_id" and (k not in id_to_row[sid] or id_to_row[sid].get(k) == "N/A"):
                                id_to_row[sid][k] = v
        except Exception as exc:
            logger.warning("Supplementary enrichment failed for %s: %s", accession, exc)

    return build_sample_table(rows, source="GEO")


def search_geo(omic, organism="human", disease=None, tissue=None, max_results=50, date_from=None, date_to=None):
    query = build_esearch_query(omic, organism, disease, tissue, date_from, date_to)
    logger.info(f"GEO query: {query}")
    api_key = os.environ.get("NCBI_API_KEY", "")
    base_params = {"api_key": api_key} if api_key else {}
    search_params = {**base_params, "db": "gds", "term": query, "retmax": max_results, "usehistory": "n"}
    resp = fetch_with_retry(ESEARCH_URL, params=search_params)
    root = ET.fromstring(resp.text)
    id_list = root.find("IdList")
    if id_list is None or not id_list.findall("Id"):
        logger.info("No GEO results found.")
        return []
    ids = [id_elem.text for id_elem in id_list.findall("Id") if id_elem.text]
    logger.info(f"Found {len(ids)} GEO IDs")
    summary_params = {**base_params, "db": "gds", "id": ",".join(ids)}
    resp = fetch_with_retry(ESUMMARY_URL, params=summary_params)
    root = ET.fromstring(resp.text)
    results = []
    for doc in root.findall("DocSum"):
        accession = _get_item_text(doc, "Accession")
        if accession and accession.startswith("GSE"):
            record = parse_geo_record(doc, omic)
            sample_table = fetch_geo_samples(accession, paper=record.get("paper"))
            if sample_table:
                record["sample_table"] = sample_table
            results.append(record)
    logger.info(f"Parsed {len(results)} GEO series records")
    return results


def main():
    parser = argparse.ArgumentParser(description="Search NCBI GEO for omics datasets")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--disease", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    parser.add_argument("--date-from", default=None)
    parser.add_argument("--date-to", default=None)
    args = parser.parse_args()
    results = search_geo(omic=args.omic, organism=args.organism, disease=args.disease, tissue=args.tissue, max_results=args.max_results, date_from=args.date_from, date_to=args.date_to)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

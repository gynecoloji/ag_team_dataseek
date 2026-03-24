#!/usr/bin/env python3
"""Search NCBI GEO for omics datasets via Entrez E-utilities."""

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET

from scripts.utils import fetch_with_retry, resolve_doi, setup_logger

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

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
            results.append(parse_geo_record(doc, omic))
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
